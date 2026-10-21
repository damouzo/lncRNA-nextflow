#!/usr/bin/env python3
"""F6.1 / F6.2: reviewer validations from the catalog + pooled expression.

F6.1 recovery_reference.tsv
    How many reference lncRNAs (catalog `origin == annotated_lncrna`) are
    recovered as expressed in this dataset (TPM >= min_tpm in >= min_samples
    samples, pooled across norm_groups). Long table: metric, value.

F6.2 tier1_vs_annotated.tsv
    Descriptive comparison of the novel Tier 1 set against reference lncRNAs:
    gene count, expressed count, fraction expressed and mean expression per set.

Expression is pooled across the per-group TPM tables because a "expressed in
>= N samples" recovery is a cohort-level descriptive claim, not a per-group DE
quantity (rule 6: no DE values are mined). Stdlib only (bin/f1_common for I/O).
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import f1_common as C  # noqa: E402

SETS = ("novel_tier1", "novel_tier2", "annotated_lncrna")


def _write(outdir, name, header, rows):
    os.makedirs(outdir, exist_ok=True)
    C.write_tsv(rows, os.path.join(outdir, name), header)


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _truthy(v):
    return str(v).strip().upper() in ("TRUE", "1", "T", "YES")


def load_expression(paths):
    """{gene_id: [tpm, ...]} pooled over every sample of every TPM table."""
    expr = {}
    for p in paths:
        header, rows = C.read_tsv(p)
        cols = [c for c in header if c != "gene_id"]
        for r in rows:
            vals = [v for v in (_num(r.get(c)) for c in cols) if v is not None]
            expr.setdefault(r["gene_id"], []).extend(vals)
    return expr


def is_expressed(vals, min_tpm, min_samples):
    return sum(1 for v in vals if v >= min_tpm) >= min_samples


def _mean(vals):
    return sum(vals) / len(vals) if vals else None


def build_recovery(genes, expr, min_tpm, min_samples):
    ref = [g for g in genes if g.get("origin") == "annotated_lncrna"]
    n_ref = len(ref)
    n_exp = sum(1 for g in ref if is_expressed(expr.get(g["gene_id"], []), min_tpm, min_samples))
    frac = ("%.4f" % (n_exp / n_ref)) if n_ref else "NA"
    n_samples = max((len(v) for v in expr.values()), default=0)
    return [
        ["min_tpm", str(min_tpm)],
        ["min_samples", str(min_samples)],
        ["n_samples_pooled", str(n_samples)],
        ["reference_lncrna_genes", str(n_ref)],
        ["reference_expressed", str(n_exp)],
        ["reference_fraction_recovered", frac],
    ]


def build_comparison(genes, expr, min_tpm, min_samples):
    buckets = {s: [] for s in SETS}
    for g in genes:
        origin, tid = g.get("origin"), g["gene_id"]
        if origin == "novel":
            key = "novel_tier1" if _truthy(g.get("tier_final")) else "novel_tier2"
        elif origin == "annotated_lncrna":
            key = "annotated_lncrna"
        else:
            continue
        buckets[key].append(tid)
    rows = []
    for s in SETS:
        ids = buckets[s]
        n = len(ids)
        n_exp = sum(1 for gid in ids if is_expressed(expr.get(gid, []), min_tpm, min_samples))
        frac = ("%.4f" % (n_exp / n)) if n else "NA"
        gene_means = [m for m in (_mean(expr.get(gid, [])) for gid in ids) if m is not None]
        mean_expr = ("%.4f" % _mean(gene_means)) if gene_means else "NA"
        rows.append([s, n, n_exp, frac, mean_expr])
    return rows


def main():
    ap = argparse.ArgumentParser(description="F6.1/F6.2 reviewer validations")
    ap.add_argument("--catalog", required=True, help="deliverables catalog_genes.tsv")
    ap.add_argument("--tpm", action="append", default=[], help="per-group <group>_tpm.tsv")
    ap.add_argument("--min-tpm", type=float, default=1.0)
    ap.add_argument("--min-samples", type=int, default=3)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    _, genes = C.read_tsv(args.catalog)
    if not genes:
        raise SystemExit("%s: empty catalog" % args.catalog)
    expr = load_expression(args.tpm) if args.tpm else {}

    _write(args.outdir, "recovery_reference.tsv", ["metric", "value"],
           build_recovery(genes, expr, args.min_tpm, args.min_samples))
    _write(args.outdir, "tier1_vs_annotated.tsv",
           ["set", "n_genes", "n_expressed", "fraction_expressed", "mean_tpm"],
           build_comparison(genes, expr, args.min_tpm, args.min_samples))
    print("recovery_analysis: wrote recovery_reference.tsv + tier1_vs_annotated.tsv")


if __name__ == "__main__":
    main()
