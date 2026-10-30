#!/usr/bin/env python3
"""Deliverables reorg tests: per-contrast enrichment split + contrasts overview.

Stdlib only. Run: python3 bin/tests/test_contrast_reorg.py
"""

import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import split_enrichment as SE            # noqa: E402
import build_contrasts_overview as BO    # noqa: E402

FAILS = []


def check(cond, msg):
    print(("  ok  " if cond else "  FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def read_rows(path):
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip() != ""]
    return lines


def test_split_enrichment():
    print("split_enrichment")
    d = tempfile.mkdtemp()
    ora = os.path.join(d, "g_enrichment_ora.tsv")
    write(ora,
          "ID\tDescription\tGeneRatio\tBgRatio\tpvalue\tp.adjust\tqvalue\tgeneID\tCount\tcontrast\n"
          "GO:1\tterm1\t1/10\t5/100\t0.01\t0.02\t0.02\tg1\t1\tA_vs_B\n"
          "GO:2\tterm2\t2/10\t6/100\t0.01\t0.03\t0.03\tg2\t2\tA_vs_B\n"
          "GO:3\tterm3\t3/10\t7/100\t0.01\t0.04\t0.04\tg3\t3\tA_vs_C\n")
    gsea = os.path.join(d, "g_enrichment_gsea.tsv")
    # header only: FUNCTIONAL_ENRICHMENT always writes the clusterProfiler header
    write(gsea,
          "ID\tDescription\tsetSize\tenrichmentScore\tNES\tpvalue\tp.adjust\t"
          "qvalue\trank\tleading_edge\tcore_enrichment\tcontrast\n")

    out = os.path.join(d, "out")
    n_ora = SE.split_table(ora, "A_vs_B", os.path.join(out, "enrichment_ora.tsv"))
    n_gsea = SE.split_table(gsea, "A_vs_B", os.path.join(out, "enrichment_gsea.tsv"))
    check(n_ora == 2, "A_vs_B keeps its 2 ORA rows")
    check(n_gsea == 0, "empty GSEA input yields 0 rows")

    ora_lines = read_rows(os.path.join(out, "enrichment_ora.tsv"))
    check(ora_lines[0].split("\t")[:3] == ["ID", "Description", "GeneRatio"],
          "ORA header preserved from the input")
    check(all("A_vs_B" in ln for ln in ora_lines[1:]) and len(ora_lines) == 3,
          "only the requested contrast rows are written")

    gsea_lines = read_rows(os.path.join(out, "enrichment_gsea.tsv"))
    check(gsea_lines[0].startswith("ID\tDescription\tsetSize"),
          "header-only GSEA input keeps its header")

    # a truly empty/missing input is an error, not a silently invented header
    try:
        SE.split_table(os.path.join(d, "nope.tsv"), "A_vs_B",
                       os.path.join(out, "missing_ora.tsv"))
        check(False, "missing input must fail")
    except SystemExit:
        check(True, "missing input fails")

    # repartition does not lose rows: per-contrast counts add up to the input
    total = n_ora + SE.split_table(ora, "A_vs_C", os.path.join(out, "c_ora.tsv"))
    check(total == 3, "per-contrast ORA rows add up to the input row count")


def test_contrasts_overview():
    print("build_contrasts_overview")
    d = tempfile.mkdtemp()
    write(os.path.join(d, "catalog_genes.tsv"),
          "gene_id\torigin\ttier_final\n"
          "N1\tnovel\tTRUE\nN2\tnovel\tFALSE\nA1\tannotated_lncrna\tTRUE\n")
    write(os.path.join(d, "staging", "c1__summary.tsv"),
          "contrast\tn_universe\tn_tested\tn_significant\tn_lncrna_full\tn_lncrna_significant\n"
          "c1\t100\t80\t5\t300\t12\n")
    write(os.path.join(d, "staging", "c1__significant.tsv"),
          "gene_id\tlog2FoldChange\tpadj\nN1\t2.0\t0.001\nA1\t1.5\t0.01\n")
    write(os.path.join(d, "staging", "c2__summary.tsv"),
          "contrast\tn_universe\tn_tested\tn_significant\nc2\t100\t70\t0\n")
    write(os.path.join(d, "staging", "c2__significant.tsv"),
          "gene_id\tlog2FoldChange\tpadj\n")

    rows = BO.build(os.path.join(d, "staging"), os.path.join(d, "catalog_genes.tsv"))
    check([r["contrast"] for r in rows] == ["c1", "c2"], "one row per contrast, sorted")
    by = {r["contrast"]: r for r in rows}
    check(by["c1"]["n_novel_tier1_significant"] == 1,
          "counts only novel Tier-1 significant genes")
    check(by["c2"]["n_novel_tier1_significant"] == 0, "0 significant -> 0 novel tier1")
    check(by["c1"]["n_lncrna_full"] == "300" and by["c1"]["n_lncrna_significant"] == "12",
          "lncRNA counts pass through from summary.tsv")
    check(by["c1"]["report"] == "contrasts/c1/c1_report.html", "report path format")

    # group-prefixed id: the report keeps the bare contrast name, no doubled prefix
    write(os.path.join(d, "staging", "grp__c3__summary.tsv"),
          "contrast\tn_universe\tn_tested\tn_significant\ngrp__c3\t100\t70\t0\n")
    write(os.path.join(d, "staging", "grp__c3__significant.tsv"),
          "gene_id\tlog2FoldChange\tpadj\n")
    by2 = {r["contrast"]: r for r in
           BO.build(os.path.join(d, "staging"), os.path.join(d, "catalog_genes.tsv"))}
    check(by2["grp__c3"]["report"] == "contrasts/grp__c3/c3_report.html",
          "group-prefixed report path drops the group prefix")


def main():
    for fn in (test_split_enrichment, test_contrasts_overview):
        fn()
    print()
    if FAILS:
        print("FAILED: %d" % len(FAILS))
        for f in FAILS:
            print("  - %s" % f)
        sys.exit(1)
    print("ALL CONTRAST REORG TESTS PASSED")


if __name__ == "__main__":
    main()
