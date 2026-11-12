#!/usr/bin/env python3
"""F1 tier_final computation (Phase A, terminal).

Applies the Tier 1 rules at the TRANSCRIPT level and aggregates to gene:

  tier_static(tx) = class tier (1a/1b) AND (multi-exonic OR tier1_allow_single_exon)
  tier_final(tx)  = tier_static(tx)          # D1: F1.4 does NOT demote

The mono-exonic rule is categorical (no derived threshold): by default a
mono-exonic transcript of any class is demoted to Tier 2 and marked
`single_exon_low_confidence`. Gene tier_static = n_tier1_transcripts >= 1;
gene single_exon_low_confidence = any demoted transcript.

F1.3 (flank) and F1.4 (antisense) are INFORMATIONAL: their columns are emitted for
traceability but do NOT gate tier_final. `antisense_status` is `evaluated` when the
element has at least one numeric antisense ratio, else `not_evaluable` (NA is never
"clean"). `x_host_leak_suspect` and `thresholds_status` travel with the table.

Also emits the F1.6 tables: per-class composition with and without
tier1_allow_single_exon, and the (now informational) threshold sensitivity.
"""

import argparse
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import f1_common as C  # noqa: E402


def flank_pass(values, thr):
    """None if unevaluable; else True when no norm_group falls below thr."""
    if not values:
        return None
    return min(values) >= thr


def antisense_risk(values, thr):
    """True when any evaluated norm_group is at or below the leakage floor."""
    return bool(values) and thr is not None and any(v <= thr for v in values)


def antisense_status(values):
    """`evaluated` when at least one numeric antisense ratio exists, else
    `not_evaluable` (an element with no usable region is never 'clean')."""
    return "evaluated" if values else "not_evaluable"


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
    antisense_risk, tier_static, tier_final). tier_final == tier_static (D1)."""
    nt = C.tier_of(class_code, parent_biotype, tier1a, tier1, parent_bios)
    if nt is None:
        return None, False, None, False, False, False
    single_exon = mono and not allow_single_exon
    static = not single_exon
    fpass = flank_pass(flank_vals.get("cand:%s" % tid, []), t_flank) if t_flank is not None else None
    arisk = antisense_risk(anti_vals.get("cand_antisense:%s" % tid, []), r_anti)
    return nt, single_exon, fpass, arisk, static, static


def count_tier1_genes(frozen, frozen_exons, gene_biotype, tier1, tier1a, parent_bios,
                      allow_single_exon):
    """Tier 1 genes under one mono-exonic setting (no thresholds involved)."""
    genes = set()
    for t in frozen["transcript"]:
        parent = C.ref_gene_of(t)
        pb = gene_biotype.get(parent) if parent else None
        nt = C.tier_of(t.get("class_code"), pb, tier1a, tier1, parent_bios)
        if nt is None:
            continue
        if len(frozen_exons.get(t["transcript_id"], [])) == 1 and not allow_single_exon:
            continue
        genes.add(t["gene_id"])
    return genes


def class_composition(frozen, frozen_exons, gene_biotype, tier1, tier1a, parent_bios,
                      allow_single_exon):
    """-> (Counter class_code -> n Tier1 transcripts,
           Counter frozenset(classes) -> n Tier1 genes) for one setting."""
    tx_by_class = Counter()
    gene_sets = {}
    for t in frozen["transcript"]:
        parent = C.ref_gene_of(t)
        pb = gene_biotype.get(parent) if parent else None
        nt = C.tier_of(t.get("class_code"), pb, tier1a, tier1, parent_bios)
        if nt is None:
            continue
        if len(frozen_exons.get(t["transcript_id"], [])) == 1 and not allow_single_exon:
            continue
        cc = t.get("class_code")
        tx_by_class[cc] += 1
        gene_sets.setdefault(t["gene_id"], set()).add(cc)
    return tx_by_class, Counter(frozenset(s) for s in gene_sets.values())


def _num(x):
    if x is None or x == "" or x == "NA":
        return None
    try:
        return float(x)
    except ValueError:
        return None


def _fmt(x):
    return "NA" if x is None else ("%.6g" % x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frozen-gtf", required=True)
    ap.add_argument("--reference-gtf", required=True)
    ap.add_argument("--flank-metrics", default="",
                    help="optional; F1.3 disabled leaves this empty")
    ap.add_argument("--antisense-metrics", required=True)
    ap.add_argument("--thresholds", required=True)
    ap.add_argument("--sensitivity-thresholds", default="")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tier1-classes", default="u,x,i")
    ap.add_argument("--tier1a-classes", default="u,x,i")
    ap.add_argument("--tier1-parent-biotypes", default="lncRNA")
    ap.add_argument("--annotated-lncrna-biotypes", default="lncRNA")
    ap.add_argument("--tier1-allow-single-exon", action="store_true")
    args = ap.parse_args()
    tier1 = set(C.parse_list(args.tier1_classes))
    tier1a = set(C.parse_list(args.tier1a_classes))
    parent_bios = set(C.parse_list(args.tier1_parent_biotypes))
    allow_single_exon = args.tier1_allow_single_exon

    for p in (args.frozen_gtf, args.reference_gtf, args.antisense_metrics,
              args.thresholds):
        C.ensure_file(p, "input")
    os.makedirs(args.outdir, exist_ok=True)

    _, thr_rows = C.read_tsv(args.thresholds)
    thr = {r["metric"]: r for r in thr_rows}
    t_flank = _num(thr["flank_ratio"]["value"]) if "flank_ratio" in thr else None
    r_anti = _num(thr["antisense_ratio"]["value"]) if "antisense_ratio" in thr else None
    status = ", ".join(sorted({r.get("status", "") for r in thr_rows if r.get("status")})) or ""

    flank_vals = {}
    if args.flank_metrics:
        _, flank_rows = C.read_tsv(args.flank_metrics)
        flank_vals = values_by_element(flank_rows, "flank_ratio_median")
    _, anti_rows = C.read_tsv(args.antisense_metrics)
    anti_vals = values_by_element(anti_rows, "ratio_median")

    # load_reference version-strips keys, so gene_biotype/tx2gene/gene_loc line up
    # with the version-stripped parent from C.ref_gene_of. The class-i sense/
    # antisense split needs the reference host gene: resolve it through
    # cmp_ref -> reference tx2gene (class i carries no ref_gene_id), the same route
    # F1 regions uses.
    ref = C.load_reference(args.reference_gtf, C.parse_list(args.annotated_lncrna_biotypes))
    gene_biotype = ref["gene_biotype"]
    ref_tx2gene = ref["tx2gene"]
    gene_loc = ref["gene_loc"]
    frozen = C.load_gtf(args.frozen_gtf, {"transcript", "exon"})
    frozen_exons = C.exon_map(frozen)

    tx_rows = []
    gene_static, gene_final, gene_novel_type, gene_single_exon = {}, {}, {}, {}
    gene_anti_status, gene_anti_suspect = {}, {}
    for t in frozen["transcript"]:
        tid = t["transcript_id"]
        parent = C.ref_gene_of(t)
        pb = gene_biotype.get(parent) if parent else None
        i_rel = C.intronic_relation(t.get("class_code"), t["strand"],
                                    C.host_gene_of(t, ref_tx2gene), gene_loc)
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
        avals = anti_vals.get("cand_antisense:%s" % tid, [])
        a_status = antisense_status(avals)
        if a_status == "evaluated":
            gene_anti_status[gene] = "evaluated"
        if arisk:
            gene_anti_suspect[gene] = True
        tx_rows.append((tid, gene, t.get("class_code"), pb, nt or "2",
                        "TRUE" if mono else "FALSE",
                        "TRUE" if single_exon else "FALSE",
                        _fmt(C.median(flank_vals.get("cand:%s" % tid, []))),
                        "NA" if fpass is None else ("TRUE" if fpass else "FALSE"),
                        a_status,
                        _fmt(C.median(avals)), len(avals),
                        "TRUE" if arisk else "FALSE", status,
                        "TRUE" if static else "FALSE",
                        "TRUE" if final else "FALSE",
                        i_rel or "NA"))

    C.write_tsv(tx_rows, os.path.join(args.outdir, "tier_final_transcripts.tsv"),
                ["transcript_id", "gene_id", "class_code", "parent_biotype",
                 "novel_type", "monoexonic", "single_exon_low_confidence",
                 "flank_ratio_median", "flank_pass", "antisense_status",
                 "antisense_ratio_median", "n_groups_evaluated",
                 "x_host_leak_suspect", "thresholds_status",
                 "tier_static", "tier_final", "i_relation"])

    gene_rows = []
    for gene in sorted(set(gene_static) | set(gene_final)):
        ns = gene_static.get(gene, 0)
        gene_rows.append((gene, gene_novel_type.get(gene, "2"),
                          "TRUE" if ns >= 1 else "FALSE", ns,
                          "TRUE" if gene_single_exon.get(gene) else "FALSE",
                          gene_anti_status.get(gene, "not_evaluable"),
                          "TRUE" if gene_anti_suspect.get(gene) else "FALSE",
                          "TRUE" if gene_final.get(gene, 0) >= 1 else "FALSE",
                          status))
    C.write_tsv(gene_rows, os.path.join(args.outdir, "tier_final_genes.tsv"),
                ["gene_id", "novel_type", "tier_static", "n_tier1_transcripts",
                 "single_exon_low_confidence", "antisense_status",
                 "x_host_leak_suspect", "tier_final", "thresholds_status"])

    # ---- F1.6 per-class composition, with and without tier1_allow_single_exon ----
    by_class_rows = []
    totals = {}
    for allow in (False, True):
        tag = "TRUE" if allow else "FALSE"
        tx_by_class, genes_by_set = class_composition(
            frozen, frozen_exons, gene_biotype, tier1, tier1a, parent_bios, allow)
        totals[allow] = sum(genes_by_set.values())
        for cc in sorted(tx_by_class):
            by_class_rows.append((tag, "transcript", cc, tx_by_class[cc]))
        for cs in sorted(genes_by_set, key=lambda s: (len(s), sorted(s))):
            by_class_rows.append((tag, "gene", "+".join(sorted(cs)), genes_by_set[cs]))
        by_class_rows.append((tag, "gene", "ALL", totals[allow]))
    C.write_tsv(by_class_rows, os.path.join(args.outdir, "tier_final_by_class.tsv"),
                ["tier1_allow_single_exon", "unit", "class_set", "n"])

    # ---- F1.6 sensitivity: thresholds are informational after D1 (tier_final is
    # categorical), so the counts here do not move; kept for the record. ----
    n_tier1 = sum(1 for r in gene_rows if r[7] == "TRUE")
    sens_rows = [("baseline", "derived", _fmt(t_flank), _fmt(r_anti), str(n_tier1))]
    if args.sensitivity_thresholds:
        _, sens = C.read_tsv(args.sensitivity_thresholds)
        for row in sens:
            sens_rows.append((row["metric"], row["percentile"], _fmt(t_flank),
                              _fmt(r_anti), str(n_tier1)))
    C.write_tsv(sens_rows, os.path.join(args.outdir, "tier_final_sensitivity.tsv"),
                ["changed", "percentile", "flank_ratio_threshold",
                 "antisense_ratio_threshold", "n_tier1_final_genes"])

    print("F1 tier_final: %d transcripts, %d Tier1 genes (allow_single_exon=%s), "
          "%d Tier1 genes (allow_single_exon=true), status=%s"
          % (len(tx_rows), totals[False], "false", totals[True], status or "ok"))


if __name__ == "__main__":
    main()
