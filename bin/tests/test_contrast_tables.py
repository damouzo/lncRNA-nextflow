#!/usr/bin/env python3
"""F4.2 tests: build_contrast_tables across several contrasts per group.
Stdlib only, no Nextflow. Validates that the LFC cutoff only affects the
`significant` flag (never the padj / DESeq2 null hypothesis).

Run: python3 bin/tests/test_contrast_tables.py
"""

import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.dirname(HERE)

FAILS = []


def check(cond, msg):
    print(("  ok  " if cond else "  FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def read_tsv(path):
    with open(path) as fh:
        lines = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
    header = lines[0]
    return [dict(zip(header, row)) for row in lines[1:]]


def build(group, ct, d, out_dir):
    universe = os.path.join(d, "universe_%s.tsv" % group)
    de = os.path.join(d, "de_%s_%s.tsv" % (group, ct))
    catalog = os.path.join(d, "catalog.tsv")
    return subprocess.run([
        sys.executable, os.path.join(BIN, "build_contrast_tables.py"),
        "--universe", universe, "--de", de, "--catalog", catalog,
        "--contrast", "%s__%s" % (group, ct), "--fdr", "0.05", "--lfc", "1.0",
        "--out", os.path.join(out_dir, "de_full.tsv"),
        "--significant-out", os.path.join(out_dir, "de_significant.tsv"),
        "--summary-out", os.path.join(out_dir, "summary.tsv"),
    ], capture_output=True, text=True)


def main():
    d = tempfile.mkdtemp()
    write(os.path.join(d, "catalog.tsv"),
          "gene_id\torigin\tclass_code\ttier_static\ttier_final\tsingle_exon_low_confidence\tantisense_status\n"
          "G1\tnovel\tu\tTRUE\tTRUE\tFALSE\tnot_evaluated\n"
          "G2\tnovel\tu\tFALSE\tFALSE\tFALSE\tnot_evaluated\n"
          "G3\tannotated_lncrna\tNA\tNA\tFALSE\tNA\tnot_evaluated\n")

    # universe = pre-filter universe (all imported genes), one file per group
    for group in ("MNC", "MSC"):
        write(os.path.join(d, "universe_%s.tsv" % group),
              "gene_id\tkeep\tfilter_status\tfilter_reason\n"
              "G1\tTRUE\tkept\t\nG2\tTRUE\tkept\t\nG3\tFALSE\tfiltered\tlow_expression\n")

    # DE tables: G1 padj<fdr and |lfc|>=1 -> significant.
    # G2 padj<fdr but |lfc|=0.3 (<lfc) -> NOT significant (cut only on the flag).
    # G3 absent from DE -> tested=FALSE, filtered_low_expression.
    for group in ("MNC", "MSC"):
        for ct in ("A_vs_B", "A_vs_C"):
            write(os.path.join(d, "de_%s_%s.tsv" % (group, ct)),
                  "gene_id\tbaseMean\tlog2FoldChange\tlfcSE\tpvalue\tpadj\n"
                  "G1\t100\t2.5\t0.3\t0.0001\t0.001\n"
                  "G2\t80\t0.3\t0.1\t0.0001\t0.001\n")

    n_runs = 0
    for group in ("MNC", "MSC"):
        for ct in ("A_vs_B", "A_vs_C"):
            out_dir = os.path.join(d, "contrasts", "%s__%s" % (group, ct))
            p = build(group, ct, d, out_dir)
            n_runs += 1
            check(p.returncode == 0, "build %s__%s runs" % (group, ct))

            full = read_tsv(os.path.join(out_dir, "de_full.tsv"))
            sig = read_tsv(os.path.join(out_dir, "de_significant.tsv"))
            summ = read_tsv(os.path.join(out_dir, "summary.tsv"))[0]

            check(len(full) == 3, "%s__%s: de_full has one row per universe gene"
                  % (group, ct))
            check(all(r["contrast"] == "%s__%s" % (group, ct) for r in full),
                  "%s__%s: contrast id is <grupo>__<contraste>" % (group, ct))
            by = {r["gene_id"]: r for r in full}
            check(by["G1"]["significant"] == "TRUE", "%s__%s: G1 significant" % (group, ct))
            check(by["G2"]["significant"] == "FALSE",
                  "%s__%s: G2 padj<fdr but |lfc|<cut -> not significant (LFC on the flag only)"
                  % (group, ct))
            check(by["G3"]["tested"] == "FALSE" and by["G3"]["not_tested_reason"] == "filtered_low_expression",
                  "%s__%s: gene absent from DE is kept as tested=FALSE" % (group, ct))
            check(int(summ["n_significant"]) == len(sig) == 1,
                  "%s__%s: summary n_significant matches de_significant" % (group, ct))
            check(by["G1"]["antisense_status"] == "not_evaluated",
                  "%s__%s: catalog annotation is joined in" % (group, ct))

    check(n_runs == 4, "two groups x two contrasts were built")

    print()
    if FAILS:
        print("FAILED: %d" % len(FAILS))
        for f in FAILS:
            print("  - %s" % f)
        sys.exit(1)
    print("ALL CONTRAST TABLE TESTS PASSED")


if __name__ == "__main__":
    main()
