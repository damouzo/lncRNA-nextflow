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
      --significant-out contrasts/<id>/de_significant.tsv \
      --summary-out contrasts/<id>/summary.tsv

The LFC cutoff is applied only to the `significant` flag here (|log2FC| >= lfc);
it must not enter DESeq2 results() as lfcThreshold, which would change the null
hypothesis and the padj values.
"""

import argparse
import csv
import os
import sys

IS_NA = {"", "NA", "NaN", "nan", "None", "NA_character_"}
ANNOT_COLS = ["origin", "class_code", "tier_static", "tier_final",
              "single_exon_low_confidence", "antisense_status"]
DE_COLS = ["baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"]


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


def build_contrast(universe_ids, catalog_by_id, de_by_id, contrast, fdr, lfc_threshold):
    rows = []
    for gid in universe_ids:
        ann = catalog_by_id.get(gid, {})
        d = de_by_id.get(gid)
        tested = d is not None
        if tested and len(d) >= 0:
            baseMean, lfc, lfcse = d.get("baseMean", ""), d.get("log2FoldChange", ""), d.get("lfcSE", "")
            pval, padj = d.get("pvalue", ""), d.get("padj", "")
            na_reason = padj_na_reason(baseMean, pval, padj)
            sig = is_significant(padj, lfc, fdr, lfc_threshold)
            reason = ""
        else:
            baseMean = lfc = lfcse = pval = padj = ""
            na_reason = ""
            sig = False
            reason = "filtered_low_expression"
        rows.append({
            "gene_id": gid, "contrast": contrast,
            "baseMean": baseMean, "log2FoldChange": lfc, "lfcSE": lfcse,
            "pvalue": pval, "padj": padj,
            "significant": "TRUE" if sig else "FALSE",
            "tested": "TRUE" if tested else "FALSE",
            "not_tested_reason": reason, "padj_na_reason": na_reason,
            **{c: ann.get(c, "") for c in ANNOT_COLS},
        })
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
    ap.add_argument("--out", required=True)
    ap.add_argument("--significant-out")
    ap.add_argument("--summary-out")
    args = ap.parse_args()

    universe = gene_ids(args.universe)
    catalog = _index_unique(args.catalog, "gene_id")
    de = _index_unique(args.de, "gene_id")

    rows = build_contrast(universe, catalog, de, args.contrast, args.fdr, args.lfc)

    cols = (["gene_id", "contrast", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj",
             "significant", "tested", "not_tested_reason", "padj_na_reason"] + ANNOT_COLS)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    n_tested = sum(1 for r in rows if r["tested"] == "TRUE")
    n_sig = sum(1 for r in rows if r["significant"] == "TRUE")
    if args.significant_out:
        sig_rows = [r for r in rows if r["significant"] == "TRUE"]
        os.makedirs(os.path.dirname(os.path.abspath(args.significant_out)), exist_ok=True)
        with open(args.significant_out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            w.writeheader()
            w.writerows(sig_rows)
    if args.summary_out:
        os.makedirs(os.path.dirname(os.path.abspath(args.summary_out)), exist_ok=True)
        with open(args.summary_out, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["contrast", "n_universe", "n_tested", "n_significant"])
            w.writerow([args.contrast, len(rows), n_tested, n_sig])

    print("build_contrast_tables [%s]: %d genes (%d tested, %d significant)"
          % (args.contrast, len(rows), n_tested, n_sig))


if __name__ == "__main__":
    main()
