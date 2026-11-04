#!/usr/bin/env python3
"""F4.2: build the per-contrast master table (`de_full`) without losing genes.

Stdlib only (csv), no pandas. The master table has one row per gene of the
quantified universe — significant or not, tested or not — so "nothing is lost" is
enforced by construction and by `assert_nrow` after every join.

Joins:
  universe (all quantified genes)  LEFT JOIN  catalog (annotation)  LEFT JOIN  DE
A gene absent from the DE table is `tested = FALSE` with
`not_tested_reason = filtered_low_expression` (DESeq2 independent filtering);
a tested gene with `padj = NA` carries `padj_na_reason`.

CLI:
  build_contrast_tables.py --universe genes.tsv --de DE.tsv --catalog catalog_genes.tsv \
      --contrast <id> --fdr 0.05 --lfc 0.0 --out contrasts/<id>/de_full.tsv \
      --significant-out contrasts/<id>/de_full_significant.tsv \
      --lncrna-out contrasts/<id>/de_lncrna_full.tsv \
      --lncrna-significant-out contrasts/<id>/de_lncrna_significant.tsv \
      --summary-out contrasts/<id>/summary.tsv

`significant_padj` (padj < FDR) is the default "significant" everywhere; the
`significant_padj_lfc1` column adds the |log2FC| >= lfc cutoff and is only used
where a table/panel names it explicitly. The `*_significant` tables ship the
`significant_padj` set; the lfc1 subset is a column inside them. The LFC cutoff is
applied only to that column here; it must not enter DESeq2 results() as
lfcThreshold, which would change the null hypothesis and the padj values.

`summary.tsv` is a one-row index of the row count of every table shipped for the
contrast: `n_universe` (de_full), `n_tested`, `n_significant` (de_full_significant),
`n_lncrna_full` (de_lncrna_full) and `n_lncrna_significant` (de_lncrna_significant).
"""

import argparse
import csv
import os
import sys

IS_NA = {"", "NA", "NaN", "nan", "None", "NA_character_"}
ANNOT_COLS = ["origin", "class_code", "tier", "tier_static", "tier_final",
              "single_exon_low_confidence", "antisense_status",
              "gene_name", "gene_biotype", "catalog_id", "locus", "lnc_position_class",
              "known_catalog_hit", "catalog_sources", "novelty_class", "n_tx_hit", "n_tx",
              "in_lncipedia", "lncipedia_best_id", "lncipedia_best_pct",
              "in_noncode", "noncode_best_id", "noncode_best_pct",
              "in_rnacentral", "rnacentral_best_id", "rnacentral_best_pct",
              "in_fantomcat", "fantomcat_best_id", "fantomcat_best_pct"]
# 4b.1: genomic position (from the catalog) and derived grouping/QC columns.
POS_COLS = ["chr", "start", "end", "strand"]
DERIVED_COLS = ["gene_id_stable", "gene_group", "contig_type", "baseMean_extreme", "genome_build"]
DE_COLS = ["baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj",
           "stat", "lfc_unshrunken", "shrink_type"]

# Canonical GRCh38 chromosomes; anything else (patches, scaffolds, alt contigs
# like KI270733.1) is non-canonical and flagged so wet lab does not chase a
# locus that may be a repeat/rRNA artefact.
CANONICAL_CONTIGS = set([str(i) for i in range(1, 23)]) | {"X", "Y", "MT", "M"}


def _is_na(v):
    return v is None or str(v) in IS_NA


def gene_group(ann):
    """4b.4: coarse group derived from biotype, not just origin (which lumps
    protein-coding, pseudogenes, TEC and TR genes under 'reference_other')."""
    origin = ann.get("origin", "")
    if origin == "novel":
        return "novel_lncRNA"
    if origin == "annotated_lncrna":
        return "annotated_lncRNA"
    if (ann.get("gene_biotype") or "").strip() == "protein_coding":
        return "protein_coding"
    return "other"


def contig_type(chrom):
    if _is_na(chrom):
        return "NA"
    return "canonical" if str(chrom).strip() in CANONICAL_CONTIGS else "noncanonical"



def read_tsv(path):
    with open(path) as fh:
        rdr = csv.reader(fh, delimiter="\t")
        header = next(rdr, None)
        if header is None:
            return [], []
        return header, [dict(zip(header, row)) for row in rdr if row]


def _num(x):
    if x is None or x in IS_NA:
        return None
    try:
        return float(x)
    except ValueError:
        return None


def gene_ids(path):
    """First column named gene_id, else the first column; order preserved, dups checked."""
    header, rows = read_tsv(path)
    key = "gene_id" if "gene_id" in header else header[0]
    ids = [r[key] for r in rows]
    if len(set(ids)) != len(ids):
        raise SystemExit("universe has duplicate gene_id: %s" % path)
    return ids


def assert_nrow(actual, expected, what):
    if actual != expected:
        raise SystemExit("assert_nrow failed after %s: %d != %d" % (what, actual, expected))


def padj_na_reason(baseMean, pvalue, padj):
    """Why a tested gene has padj = NA, from the fields we have."""
    if padj is not None and str(padj) not in IS_NA:
        return ""
    bm = _num(baseMean)
    if bm is not None and bm == 0:
        return "baseMean_zero"
    return "independent_filtering"


def is_significant(padj, log2fc, fdr, lfc_threshold):
    p = _num(padj)
    f = _num(log2fc)
    if p is None or f is None:
        return False
    return p < fdr and abs(f) >= lfc_threshold


def is_significant_padj(padj, fdr):
    """4a.3: padj-only significance (the definition the report figures used)."""
    p = _num(padj)
    return p is not None and p < fdr


def build_contrast(universe_ids, catalog_by_id, de_by_id, contrast, fdr, lfc_threshold,
                   extreme_basemean, genome_build):
    rows = []
    for gid in universe_ids:
        ann = catalog_by_id.get(gid, {})
        d = de_by_id.get(gid)
        tested = d is not None
        if tested and len(d) >= 0:
            baseMean, lfc, lfcse = d.get("baseMean", ""), d.get("log2FoldChange", ""), d.get("lfcSE", "")
            pval, padj = d.get("pvalue", ""), d.get("padj", "")
            stat = d.get("stat", "")
            lfc_unshrunken = d.get("lfc_unshrunken", "")
            shrink_type = d.get("shrink_type", "")
            na_reason = padj_na_reason(baseMean, pval, padj)
            sig = is_significant(padj, lfc, fdr, lfc_threshold)
            sig_padj = is_significant_padj(padj, fdr)
            reason = ""
        else:
            baseMean = lfc = lfcse = pval = padj = ""
            stat = lfc_unshrunken = shrink_type = ""
            na_reason = ""
            sig = False
            sig_padj = False
            reason = "filtered_low_expression"
        bm = _num(baseMean)
        row = {
            "gene_id": gid, "contrast": contrast,
            "baseMean": baseMean, "log2FoldChange": lfc, "lfcSE": lfcse,
            "pvalue": pval, "padj": padj,
            "stat": stat, "lfc_unshrunken": lfc_unshrunken, "shrink_type": shrink_type,
            "significant_padj": "TRUE" if sig_padj else "FALSE",
            "significant_padj_lfc1": "TRUE" if sig else "FALSE",
            "tested": "TRUE" if tested else "FALSE",
            "not_tested_reason": reason, "padj_na_reason": na_reason,
            **{c: ann.get(c, "") for c in ANNOT_COLS},
        }
        # 4b.1 / 4a.4: positions from the catalog + derived grouping and QC flags.
        row["chr"] = ann.get("chromosome", "")
        row["start"] = ann.get("start", "")
        row["end"] = ann.get("end", "")
        row["strand"] = ann.get("strand", "")
        # 4b.1 (decision b): one stable id column for every gene — the catalog id
        # for novel genes, the reference gene_id otherwise.
        cid = ann.get("catalog_id")
        row["gene_id_stable"] = cid if not _is_na(cid) else gid
        row["gene_group"] = gene_group(ann)
        row["contig_type"] = contig_type(ann.get("chromosome"))
        row["baseMean_extreme"] = (
            "TRUE" if (bm is not None and extreme_basemean and bm >= extreme_basemean) else "FALSE")
        row["genome_build"] = genome_build
        rows.append(row)
    assert_nrow(len(rows), len(universe_ids), "universe x catalog x DE join")
    return rows


def _index_unique(path, key):
    header, rows = read_tsv(path)
    out = {}
    for r in rows:
        k = r.get(key)
        if k in out:
            raise SystemExit("duplicate %s '%s' in %s" % (key, k, path))
        out[k] = r
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--universe", required=True, help="TSV with a gene_id column (all quantified genes)")
    ap.add_argument("--de", required=True, help="DESeq2 results TSV (tested genes)")
    ap.add_argument("--catalog", required=True, help="catalog gene TSV with annotation columns")
    ap.add_argument("--contrast", required=True)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--lfc", type=float, default=0.0)
    ap.add_argument("--extreme-basemean", type=float, default=100000.0,
                    help="4a.4: baseMean at/above which baseMean_extreme is TRUE")
    ap.add_argument("--genome-build", default="",
                    help="4b.1: genome build recorded on every row (e.g. GRCh38)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--significant-out")
    ap.add_argument("--lncrna-out",
                    help="4b.4: unfiltered lncRNA-only table (annotated + novel)")
    ap.add_argument("--lncrna-significant-out",
                    help="4b.5: lncRNAs with significant_padj == TRUE (padj < FDR)")
    ap.add_argument("--summary-out")
    args = ap.parse_args()

    universe = gene_ids(args.universe)
    catalog = _index_unique(args.catalog, "gene_id")
    de = _index_unique(args.de, "gene_id")

    rows = build_contrast(universe, catalog, de, args.contrast, args.fdr, args.lfc,
                          args.extreme_basemean, args.genome_build)

    cols = (["gene_id", "contrast", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj",
             "stat", "lfc_unshrunken", "shrink_type",
             "significant_padj", "significant_padj_lfc1",
             "tested", "not_tested_reason", "padj_na_reason"] + ANNOT_COLS
            + POS_COLS + DERIVED_COLS)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    n_tested = sum(1 for r in rows if r["tested"] == "TRUE")
    # "significant" = padj < FDR (significant_padj); the lfc1 subset is a column.
    n_sig = sum(1 for r in rows if r["significant_padj"] == "TRUE")
    n_lnc = sum(1 for r in rows if r["gene_group"] in ("annotated_lncRNA", "novel_lncRNA"))
    n_lnc_sig = sum(1 for r in rows
                    if r["gene_group"] in ("annotated_lncRNA", "novel_lncRNA")
                    and r["significant_padj"] == "TRUE")
    if args.significant_out:
        sig_rows = [r for r in rows if r["significant_padj"] == "TRUE"]
        os.makedirs(os.path.dirname(os.path.abspath(args.significant_out)), exist_ok=True)
        with open(args.significant_out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            w.writeheader()
            w.writerows(sig_rows)
    if args.lncrna_out:
        # 4b.4: every lncRNA (annotated + novel), significant or not, so the wet
        # lab gets a complete table to open in Excel without the coding "background".
        lnc_rows = [r for r in rows if r["gene_group"] in ("annotated_lncRNA", "novel_lncRNA")]
        os.makedirs(os.path.dirname(os.path.abspath(args.lncrna_out)), exist_ok=True)
        with open(args.lncrna_out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            w.writeheader()
            w.writerows(lnc_rows)
    if args.lncrna_significant_out:
        # 4b.5: significant lncRNAs (annotated + novel), significant = significant_padj.
        lnc_sig = [r for r in rows
                   if r["gene_group"] in ("annotated_lncRNA", "novel_lncRNA")
                   and r["significant_padj"] == "TRUE"]
        os.makedirs(os.path.dirname(os.path.abspath(args.lncrna_significant_out)), exist_ok=True)
        with open(args.lncrna_significant_out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            w.writeheader()
            w.writerows(lnc_sig)
    if args.summary_out:
        # One row indexing the row counts of every table shipped for this contrast:
        # n_universe -> de_full, n_significant -> de_full_significant,
        # n_lncrna_full -> de_lncrna_full, n_lncrna_significant -> de_lncrna_significant.
        os.makedirs(os.path.dirname(os.path.abspath(args.summary_out)), exist_ok=True)
        with open(args.summary_out, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["contrast", "n_universe", "n_tested", "n_significant",
                        "n_lncrna_full", "n_lncrna_significant"])
            w.writerow([args.contrast, len(rows), n_tested, n_sig, n_lnc, n_lnc_sig])

    print("build_contrast_tables [%s]: %d genes (%d tested, %d significant)"
          % (args.contrast, len(rows), n_tested, n_sig))


if __name__ == "__main__":
    main()
