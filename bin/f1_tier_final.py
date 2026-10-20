#!/usr/bin/env python3
"""F1 tier_final computation (Phase A, terminal).

Applies the Tier 1 rules at the TRANSCRIPT level and aggregates to gene:

  tier_static(tx) = class tier (1a/1b) AND (multi-exonic OR tier1_allow_single_exon)
  tier_final(tx)  = tier_static(tx) AND not antisense-leakage flag

The mono-exonic rule is categorical (no derived threshold): by default a
mono-exonic transcript of any class is demoted to Tier 2 and marked
`single_exon_low_confidence`. Gene tier_static = n_tier1_transcripts >= 1;
gene single_exon_low_confidence = any demoted transcript.

F1.3 (flank coverage) is exploratory and NOT applied to tier_final: its columns
are reported when metrics are present, otherwise NA. tier_final depends only on
the class/structure rule and the F1.4 antisense flag.

Also emits the F1.6 sensitivity table (alternative antisense percentiles, plus
flank percentiles when flank metrics are present).
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import f1_common as C  # noqa: E402


def flank_pass(values, thr):
    """None if unevaluable; else True when no norm_group falls below thr."""
    if not values:
        return None
    return min(values) >= thr


def antisense_risk(values, thr):
    return bool(values) and any(v <= thr for v in values)


def values_by_element(rows, value_col):
    out = {}
    for r in rows:
        v = _num(r.get(value_col))
        if v is not None:
            out.setdefault(r["element_id"], []).append(v)
    return out


def evaluate_tx(tid, class_code, parent_biotype, mono, flank_vals, anti_vals,
                t_flank, r_anti, tier1, tier1a, parent_bios, allow_single_exon):
    """-> (novel_type|None, single_exon_low_confidence, flank_pass|None,
    antisense_risk, tier_static, tier_final)."""
    nt = C.tier_of(class_code, parent_biotype, tier1a, tier1, parent_bios)
    if nt is None:
        return None, False, None, False, False, False
    single_exon = mono and not allow_single_exon
    static = not single_exon
    fpass = flank_pass(flank_vals.get("cand:%s" % tid, []), t_flank) if t_flank is not None else None
    arisk = antisense_risk(anti_vals.get("cand_antisense:%s" % tid, []), r_anti) if r_anti is not None else False
    final = static and not arisk
    return nt, single_exon, fpass, arisk, static, final


def count_final_genes(frozen, frozen_exons, gene_biotype, tier1, tier1a,
                      parent_bios, flank_vals, anti_vals, t_flank, r_anti,
                      allow_single_exon):
    genes = set()
    for t in frozen["transcript"]:
        tid = t["transcript_id"]
        parent = C.ref_gene_of(t)
        pb = gene_biotype.get(parent) if parent else None
        _, _, _, _, _, final = evaluate_tx(
            tid, t.get("class_code"), pb, len(frozen_exons.get(tid, [])) == 1,
            flank_vals, anti_vals, t_flank, r_anti, tier1, tier1a, parent_bios,
            allow_single_exon)
        if final:
            genes.add(t["gene_id"])
    return len(genes)


def _num(x):
    if x is None or x == "" or x == "NA":
        return None
    try:
        return float(x)
    except ValueError:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frozen-gtf", required=True)
    ap.add_argument("--reference-gtf", required=True)
    ap.add_argument("--flank-metrics", required=True)
    ap.add_argument("--antisense-metrics", required=True)
    ap.add_argument("--thresholds", required=True)
    ap.add_argument("--sensitivity-thresholds", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tier1-classes", default="u,x,i,j,o")
    ap.add_argument("--tier1a-classes", default="u,x,i")
    ap.add_argument("--tier1-parent-biotypes", default="lncRNA")
    ap.add_argument("--tier1-allow-single-exon", action="store_true")
    args = ap.parse_args()
    tier1 = set(C.parse_list(args.tier1_classes))
    tier1a = set(C.parse_list(args.tier1a_classes))
    parent_bios = set(C.parse_list(args.tier1_parent_biotypes))
    allow_single_exon = args.tier1_allow_single_exon

    for p in (args.frozen_gtf, args.reference_gtf, args.flank_metrics,
              args.antisense_metrics, args.thresholds, args.sensitivity_thresholds):
        C.ensure_file(p, "input")
    os.makedirs(args.outdir, exist_ok=True)

    _, thr_rows = C.read_tsv(args.thresholds)
    thr = {r["metric"]: r for r in thr_rows}
    t_flank = _num(thr["flank_ratio"]["value"]) if "flank_ratio" in thr else None
    r_anti = _num(thr["antisense_ratio"]["value"]) if "antisense_ratio" in thr else None
    status = ", ".join(sorted({r.get("status", "") for r in thr_rows if r.get("status")})) or ""

    _, flank_rows = C.read_tsv(args.flank_metrics)
    _, anti_rows = C.read_tsv(args.antisense_metrics)
    flank_vals = values_by_element(flank_rows, "flank_ratio_median")
    anti_vals = values_by_element(anti_rows, "ratio_median")

    reference = C.load_gtf(args.reference_gtf, {"gene"})
    # Strip versions so the key matches the version-stripped parent from
    # C.ref_gene_of (consistent with load_reference and the R implementation).
    gene_biotype = {C.strip_version(g["gene_id"]): g.get("gene_biotype")
                    for g in reference["gene"]}
    frozen = C.load_gtf(args.frozen_gtf, {"transcript", "exon"})
    frozen_exons = C.exon_map(frozen)

    tx_rows = []
    gene_static, gene_final, gene_novel_type, gene_single_exon = {}, {}, {}, {}
    for t in frozen["transcript"]:
        tid = t["transcript_id"]
        parent = C.ref_gene_of(t)
        pb = gene_biotype.get(parent) if parent else None
        mono = len(frozen_exons.get(tid, [])) == 1
        nt, single_exon, fpass, arisk, static, final = evaluate_tx(
            tid, t.get("class_code"), pb, mono, flank_vals, anti_vals,
            t_flank, r_anti, tier1, tier1a, parent_bios, allow_single_exon)
        gene = t["gene_id"]
        gene_static[gene] = gene_static.get(gene, 0) + (1 if static else 0)
        gene_final[gene] = gene_final.get(gene, 0) + (1 if final else 0)
        if single_exon:
            gene_single_exon[gene] = True
        if nt == "1a":
            gene_novel_type[gene] = "1a"
        elif nt == "1b" and gene_novel_type.get(gene) != "1a":
            gene_novel_type[gene] = "1b"
        tx_rows.append((tid, gene, t.get("class_code"), pb, nt or "2",
                        "TRUE" if mono else "FALSE",
                        "TRUE" if single_exon else "FALSE",
                        _fmt(C.median(flank_vals.get("cand:%s" % tid, []))),
                        "NA" if fpass is None else ("TRUE" if fpass else "FALSE"),
                        _fmt(C.median(anti_vals.get("cand_antisense:%s" % tid, []))),
                        "TRUE" if arisk else "FALSE",
                        "TRUE" if static else "FALSE",
                        "TRUE" if final else "FALSE"))

    C.write_tsv(tx_rows, os.path.join(args.outdir, "tier_final_transcripts.tsv"),
                ["transcript_id", "gene_id", "class_code", "parent_biotype",
                 "novel_type", "monoexonic", "single_exon_low_confidence",
                 "flank_ratio_median", "flank_pass", "antisense_ratio_median",
                 "antisense_risk", "tier_static", "tier_final"])

    gene_rows = []
    for gene in sorted(set(gene_static) | set(gene_final)):
        ns = gene_static.get(gene, 0)
        nf = gene_final.get(gene, 0)
        gene_rows.append((gene, gene_novel_type.get(gene, "2"),
                          "TRUE" if ns >= 1 else "FALSE", ns, nf,
                          "TRUE" if gene_single_exon.get(gene) else "FALSE",
                          "TRUE" if nf >= 1 else "FALSE", status))
    C.write_tsv(gene_rows, os.path.join(args.outdir, "tier_final_genes.tsv"),
                ["gene_id", "novel_type", "tier_static", "n_tier1_transcripts",
                 "n_tier1_transcripts_final", "single_exon_low_confidence",
                 "tier_final", "thresholds_status"])

    # ---- F1.6 sensitivity: vary one threshold at a time ----
    _, sens = C.read_tsv(args.sensitivity_thresholds)
    sens_rows = [("baseline", "derived", _fmt(t_flank), _fmt(r_anti),
                  str(sum(1 for r in gene_rows if r[6] == "TRUE")))]
    for row in sens:
        metric = row["metric"]
        alt = _num(row["value"])
        if metric == "flank_ratio":
            if t_flank is None:
                continue
            n = count_final_genes(frozen, frozen_exons, gene_biotype, tier1, tier1a,
                                  parent_bios, flank_vals, anti_vals, alt, r_anti,
                                  allow_single_exon)
            sens_rows.append((metric, row["percentile"], _fmt(alt), _fmt(r_anti), str(n)))
        else:
            n = count_final_genes(frozen, frozen_exons, gene_biotype, tier1, tier1a,
                                  parent_bios, flank_vals, anti_vals, t_flank, alt,
                                  allow_single_exon)
            sens_rows.append((metric, row["percentile"], _fmt(t_flank), _fmt(alt), str(n)))
    C.write_tsv(sens_rows, os.path.join(args.outdir, "tier_final_sensitivity.tsv"),
                ["changed", "percentile", "flank_ratio_threshold",
                 "antisense_ratio_threshold", "n_tier1_final_genes"])
    print("F1 tier_final: %d transcripts, %d Tier1 genes (status=%s)"
          % (len(tx_rows), sum(1 for r in gene_rows if r[6] == "TRUE"), status or "ok"))


def _fmt(x):
    return "NA" if x is None else ("%.6g" % x)


if __name__ == "__main__":
    main()
