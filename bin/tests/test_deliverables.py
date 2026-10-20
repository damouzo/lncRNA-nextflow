#!/usr/bin/env python3
"""F4.0 (contract validator) + F4.2 (master-table builder) tests. Stdlib only.

Run: python3 bin/tests/test_deliverables.py
"""

import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import validate_deliverables as V   # noqa: E402
import build_contrast_tables as B   # noqa: E402

FAILS = []


def check(cond, msg):
    print(("  ok  " if cond else "  FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


CONTRACT = {
    "contract": "t", "version": "1",
    "files": [
        {"path": "contrasts/<id>/de_full.tsv", "kind": "table", "min_rows": 1,
         "columns": [{"name": "gene_id", "type": "string"},
                     {"name": "padj", "type": "float"},
                     {"name": "tested", "type": "bool"}]},
        {"path": "catalog/catalog_genes.tsv", "kind": "table", "min_rows": 1,
         "columns": [{"name": "gene_id", "type": "string"},
                     {"name": "n_tier1_transcripts", "type": "int"}]},
    ],
}


def test_validator():
    print("validator")
    d = tempfile.mkdtemp()
    write(os.path.join(d, "contrasts/c1/de_full.tsv"),
          "gene_id\tpadj\ttested\nG1\t0.01\tTRUE\nG2\tNA\tFALSE\n")
    write(os.path.join(d, "catalog/catalog_genes.tsv"),
          "gene_id\tn_tier1_transcripts\nG1\t2\n")
    check(V.validate(d, CONTRACT) == [], "valid deliverables pass")

    # missing column
    d2 = tempfile.mkdtemp()
    write(os.path.join(d2, "contrasts/c1/de_full.tsv"), "gene_id\ttested\nG1\tTRUE\n")
    write(os.path.join(d2, "catalog/catalog_genes.tsv"), "gene_id\tn_tier1_transcripts\nG1\t2\n")
    errs = V.validate(d2, CONTRACT)
    check(any("missing columns" in e and "padj" in e for e in errs), "missing column reported")

    # wrong type
    d3 = tempfile.mkdtemp()
    write(os.path.join(d3, "contrasts/c1/de_full.tsv"), "gene_id\tpadj\ttested\nG1\tNA\n")
    errs = V.validate(d3, CONTRACT)
    check(any("row shorter than header" in e for e in errs), "short row reported")

    d4 = tempfile.mkdtemp()
    write(os.path.join(d4, "contrasts/c1/de_full.tsv"), "gene_id\tpadj\ttested\nG1\tabc\tTRUE\n")
    write(os.path.join(d4, "catalog/catalog_genes.tsv"), "gene_id\tn_tier1_transcripts\nG1\t2\n")
    errs = V.validate(d4, CONTRACT)
    check(any("not float" in e for e in errs), "bad float type reported")

    # missing file
    d5 = tempfile.mkdtemp()
    errs = V.validate(d5, CONTRACT)
    check(any("no file matches" in e for e in errs), "missing contract file reported")


def test_status_helpers():
    print("status helpers")
    check(B.padj_na_reason(10, 0.2, "NA") == "independent_filtering", "padj NA -> independent_filtering")
    check(B.padj_na_reason(0, "NA", "NA") == "baseMean_zero", "baseMean 0 -> baseMean_zero")
    check(B.padj_na_reason(10, 0.2, 0.03) == "", "padj present -> no reason")
    check(B.is_significant(0.01, 2.0, 0.05, 0.0) is True, "padj<fdr significant")
    check(B.is_significant(0.5, 2.0, 0.05, 0.0) is False, "padj>=fdr not significant")
    check(B.is_significant(0.01, 0.1, 0.05, 1.0) is False, "lfc below threshold not significant")


def test_builder():
    print("contrast builder")
    d = tempfile.mkdtemp()
    write(os.path.join(d, "universe.tsv"), "gene_id\nA\nB\nC\nD\n")
    write(os.path.join(d, "de.tsv"),
          "gene_id\tbaseMean\tlog2FoldChange\tlfcSE\tpvalue\tpadj\n"
          "A\t100\t2.0\t0.3\t0.001\t0.01\n"
          "B\t50\t0.1\t0.2\t0.4\t0.5\n"
          "C\t20\t\t0.2\tNA\tNA\n")
    write(os.path.join(d, "catalog_genes.tsv"),
          "gene_id\torigin\tclass_code\ttier_static\ttier_final\tsingle_exon_low_confidence\tantisense_status\n"
          "A\tnovel\tx\tTRUE\tTRUE\tFALSE\tevaluated\n"
          "B\tnovel\ti\tFALSE\tFALSE\tFALSE\tnot_evaluable\n"
          "C\tnovel\tu\tTRUE\tTRUE\tTRUE\tnot_evaluable\n"
          "D\tnovel\tu\tFALSE\tFALSE\tFALSE\tnot_evaluable\n")

    universe = B.gene_ids(os.path.join(d, "universe.tsv"))
    check(universe == ["A", "B", "C", "D"], "universe order preserved")
    catalog = B._index_unique(os.path.join(d, "catalog_genes.tsv"), "gene_id")
    de = B._index_unique(os.path.join(d, "de.tsv"), "gene_id")
    rows = B.build_contrast(universe, catalog, de, "c1", 0.05, 0.0)
    by = {r["gene_id"]: r for r in rows}
    check(len(rows) == 4, "de_full has one row per universe gene")
    check(by["A"]["significant"] == "TRUE" and by["A"]["tested"] == "TRUE", "A tested+significant")
    check(by["B"]["significant"] == "FALSE", "B tested not significant")
    check(by["C"]["padj_na_reason"] == "independent_filtering", "C tested padj NA reason")
    check(by["D"]["tested"] == "FALSE" and by["D"]["not_tested_reason"] == "filtered_low_expression",
          "D absent -> tested FALSE + reason")
    check(by["A"]["tier_final"] == "TRUE" and by["A"]["class_code"] == "x", "annotation joined")

    # duplicate universe fails
    write(os.path.join(d, "dup.tsv"), "gene_id\nA\nA\n")
    try:
        B.gene_ids(os.path.join(d, "dup.tsv"))
        check(False, "duplicate universe gene_id must fail")
    except SystemExit:
        check(True, "duplicate universe gene_id fails")

    # assert_nrow guards
    try:
        B.assert_nrow(3, 4, "test")
        check(False, "assert_nrow must fail on mismatch")
    except SystemExit:
        check(True, "assert_nrow fails on mismatch")


def test_cli_end_to_end():
    print("builder CLI")
    d = tempfile.mkdtemp()
    write(os.path.join(d, "u.tsv"), "gene_id\nA\nB\n")
    write(os.path.join(d, "de.tsv"), "gene_id\tlog2FoldChange\tpadj\nA\t1.0\t0.001\n")
    write(os.path.join(d, "cat.tsv"), "gene_id\tclass_code\ttier_static\ttier_final\nA\tx\tTRUE\tTRUE\nB;u\tFALSE\tFALSE\n".replace("B;u", "B\tu"))
    out = os.path.join(d, "contrasts/c1/de_full.tsv")
    sumout = os.path.join(d, "contrasts/c1/summary.tsv")
    old = sys.argv
    sys.argv = ["build_contrast_tables.py", "--universe", os.path.join(d, "u.tsv"),
                "--de", os.path.join(d, "de.tsv"), "--catalog", os.path.join(d, "cat.tsv"),
                "--contrast", "c1", "--out", out, "--summary-out", sumout]
    try:
        B.main()
        ok = True
    finally:
        sys.argv = old
    check(ok and os.path.exists(out) and os.path.exists(sumout), "CLI writes de_full and summary")
    _, srows = B.read_tsv(sumout)
    check(srows[0]["n_universe"] == "2" and srows[0]["n_tested"] == "1" and srows[0]["n_significant"] == "1",
          "summary counts correct")


def main():
    for fn in (test_validator, test_status_helpers, test_builder, test_cli_end_to_end):
        fn()
    print()
    if FAILS:
        print("FAILED: %d" % len(FAILS))
        for f in FAILS:
            print("  - %s" % f)
        sys.exit(1)
    print("ALL DELIVERABLES TESTS PASSED")


if __name__ == "__main__":
    main()
