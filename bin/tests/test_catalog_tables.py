#!/usr/bin/env python3
"""F4.3 / F4.3b (catalog deliverable) tests. Stdlib only.

Run: python3 bin/tests/test_catalog_tables.py
"""

import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import build_catalog_tables as B  # noqa: E402

FAILS = []


def check(cond, msg):
    print(("  ok  " if cond else "  FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def make_fixture(d):
    write(os.path.join(d, "gene_catalog.tsv"),
          "gene_id\tgene_name\tchromosome\tstart\tend\tstrand\tgene_biotype\torigin\tclass_code\t"
          "n_tier1_transcripts\ttier_static\tsingle_exon_low_confidence\tnovel_type\n"
          "G1\tNA\t1\t1\t100\t+\tNA\tnovel\tu\t2\tTRUE\tFALSE\t1a\n"
          "G2\tNAME\t1\t200\t300\t-\tlncRNA\tannotated_lncrna\tNA\tNA\tNA\tNA\tNA\n")
    write(os.path.join(d, "tier_final_genes.tsv"),
          "gene_id\tnovel_type\ttier_static\tn_tier1_transcripts\tsingle_exon_low_confidence\t"
          "antisense_status\tx_host_leak_suspect\ttier_final\tthresholds_status\n"
          "G1\t1a\tTRUE\t2\tFALSE\tevaluated\tFALSE\tTRUE\tthresholds_not_specific\n")
    write(os.path.join(d, "tier_final_transcripts.tsv"),
          "transcript_id\tgene_id\tclass_code\tparent_biotype\tnovel_type\tmonoexonic\t"
          "single_exon_low_confidence\tantisense_status\tx_host_leak_suspect\tn_groups_evaluated\t"
          "thresholds_status\ttier_static\ttier_final\n"
          "T1\tG1\tu\tNA\t1a\tFALSE\tFALSE\tevaluated\tFALSE\t1\tthresholds_not_specific\tTRUE\tTRUE\n")
    write(os.path.join(d, "gffcmp.annotated.gtf"),
          'chr1\ts\ttranscript\t1\t300\t.\t+\t.\ttranscript_id "T1"; class_code "u";\n'
          'chr1\ts\texon\t1\t300\t.\t+\t.\ttranscript_id "T1";\n'
          'chr1\ts\ttranscript\t1\t250\t.\t+\t.\ttranscript_id "T2"; class_code "u";\n'
          'chr1\ts\texon\t1\t250\t.\t+\t.\ttranscript_id "T2";\n'
          'chr1\ts\ttranscript\t1\t500\t.\t+\t.\ttranscript_id "T3"; class_code "p";\n'
          'chr1\ts\texon\t1\t500\t.\t+\t.\ttranscript_id "T3";\n'
          'chr1\ts\ttranscript\t1\t500\t.\t+\t.\ttranscript_id "T4"; class_code "u";\n'
          'chr1\ts\texon\t1\t500\t.\t+\t.\ttranscript_id "T4";\n')
    write(os.path.join(d, "recurrence_filtered.tsv"),
          "ID\tconsensus_noncoding\tpasses_recurrence\n"
          "T1\tTRUE\tTRUE\n"
          "T2\tTRUE\tFALSE\n"
          "T4\tFALSE\tTRUE\n")
    write(os.path.join(d, "frozen.gtf"),
          'chr1\ts\ttranscript\t1\t100\t.\t+\t.\ttranscript_id "T1"; gene_id "G1";\n')
    # Kept set of TRANSCRIPT_LENGTH_FILTER: class kept + length >= 200 -> T1, T2, T4
    write(os.path.join(d, "length_filtered.gtf"),
          'chr1\ts\ttranscript\t1\t300\t.\t+\t.\ttranscript_id "T1";\n'
          'chr1\ts\texon\t1\t300\t.\t+\t.\ttranscript_id "T1";\n'
          'chr1\ts\ttranscript\t1\t250\t.\t+\t.\ttranscript_id "T2";\n'
          'chr1\ts\texon\t1\t250\t.\t+\t.\ttranscript_id "T2";\n'
          'chr1\ts\ttranscript\t1\t500\t.\t+\t.\ttranscript_id "T4";\n'
          'chr1\ts\texon\t1\t500\t.\t+\t.\ttranscript_id "T4";\n')
    write(os.path.join(d, "length_filtered_short.gtf"),
          'chr1\ts\ttranscript\t1\t300\t.\t+\t.\ttranscript_id "T1";\n')


def test_build():
    d = tempfile.mkdtemp()
    make_fixture(d)
    out = os.path.join(d, "out")
    tier_final = B.index_unique(os.path.join(d, "tier_final_genes.tsv"), "gene_id")
    n_genes, n_tier1 = B.build_catalog_genes(os.path.join(d, "gene_catalog.tsv"), tier_final, out)
    check(n_genes == 2 and n_tier1 == 1, "catalog_genes: 2 genes, 1 tier_final")
    _, genes = B.read_tsv(os.path.join(out, "catalog_genes.tsv"))
    g2 = [g for g in genes if g["gene_id"] == "G2"][0]
    check(g2["tier_final"] == "FALSE" and g2["n_tier1_transcripts"] == "0",
          "reference gene defaults tier_final FALSE / n_tier1 0")
    check(g2["antisense_status"] == "not_evaluated", "reference gene antisense not_evaluated")

    n_tx = B.build_catalog_tx(os.path.join(d, "tier_final_transcripts.tsv"), out)
    check(n_tx == 1, "catalog_tx: 1 transcript")

    n_cand, rows = B.build_candidates_all(os.path.join(d, "gffcmp.annotated.gtf"),
                                          os.path.join(d, "recurrence_filtered.tsv"),
                                          os.path.join(d, "frozen.gtf"), 200, out)
    by = {r["transcript_id"]: r for r in rows}
    check(n_cand == 4, "candidates_all: 4 candidates")
    check(by["T1"]["filter_status"] == "frozen" and by["T1"]["frozen"] == "TRUE",
          "T1 frozen")
    check(by["T2"]["filter_reason"] == "low_recurrence", "T2 low_recurrence")
    check(by["T3"]["filter_reason"] == "artifact_class", "T3 artifact_class (p)")
    check(by["T4"]["filter_reason"] == "coding", "T4 coding")
    check(by["T3"]["length"] == 500, "length summed over exons")

    # Parity vs TRANSCRIPT_LENGTH_FILTER's kept set
    B.build_candidates_all(os.path.join(d, "gffcmp.annotated.gtf"),
                           os.path.join(d, "recurrence_filtered.tsv"),
                           os.path.join(d, "frozen.gtf"), 200, out,
                           os.path.join(d, "length_filtered.gtf"))
    check(True, "parity passes when replica kept set matches the length filter")

    try:
        B.build_candidates_all(os.path.join(d, "gffcmp.annotated.gtf"),
                               os.path.join(d, "recurrence_filtered.tsv"),
                               os.path.join(d, "frozen.gtf"), 200, out,
                               os.path.join(d, "length_filtered_short.gtf"))
        check(False, "parity mismatch must raise SystemExit")
    except SystemExit:
        check(True, "parity mismatch raises SystemExit")


def test_position_class_wired():
    d = tempfile.mkdtemp()
    write(os.path.join(d, "gene_catalog.tsv"),
          "gene_id\tgene_name\tchromosome\tstart\tend\tstrand\tgene_biotype\torigin\tclass_code\t"
          "n_tier1_transcripts\ttier_static\tsingle_exon_low_confidence\tnovel_type\n"
          "CG1\tNA\t1\t1000\t2000\t+\tprotein_coding\treference_other\tNA\tNA\tNA\tNA\tNA\n"
          "G1\tNA\t1\t1500\t1600\t-\tNA\tnovel\tx\t1\tTRUE\tFALSE\t1a\n")
    write(os.path.join(d, "tier_final_genes.tsv"),
          "gene_id\ttier_final\tantisense_status\n"
          "G1\tTRUE\tevaluated\n")
    out = os.path.join(d, "out")
    tier_final = B.index_unique(os.path.join(d, "tier_final_genes.tsv"), "gene_id")
    B.build_catalog_genes(os.path.join(d, "gene_catalog.tsv"), tier_final, out)
    _, genes = B.read_tsv(os.path.join(out, "catalog_genes.tsv"))
    by = {g["gene_id"]: g for g in genes}
    check(by["G1"]["lnc_position_class"] == "antisense",
          "build_catalog_genes adds lnc_position_class (antisense)")
    check(by["G1"]["lnc_position_gene_id"] == "CG1", "position gene_id recorded")
    check(by["CG1"]["lnc_position_class"] == "NA", "reference coding row gets NA")


def test_cli():
    d = tempfile.mkdtemp()
    make_fixture(d)
    out = os.path.join(d, "out")
    old = sys.argv
    sys.argv = ["build_catalog_tables.py",
                "--gene-catalog", os.path.join(d, "gene_catalog.tsv"),
                "--tier-final-genes", os.path.join(d, "tier_final_genes.tsv"),
                "--tier-final-transcripts", os.path.join(d, "tier_final_transcripts.tsv"),
                "--gffcmp-gtf", os.path.join(d, "gffcmp.annotated.gtf"),
                "--recurrence", os.path.join(d, "recurrence_filtered.tsv"),
                "--frozen-gtf", os.path.join(d, "frozen.gtf"),
                "--outdir", out]
    try:
        B.main()
        ok = True
    finally:
        sys.argv = old
    check(ok, "CLI runs")
    for f in ("catalog_genes.tsv", "catalog_tx.tsv", "candidates_all.tsv"):
        check(os.path.exists(os.path.join(out, f)), "CLI wrote %s" % f)


def main():
    for fn in (test_build, test_position_class_wired, test_cli):
        fn()
    print()
    if FAILS:
        print("FAILED: %d" % len(FAILS))
        for f in FAILS:
            print("  - %s" % f)
        sys.exit(1)
    print("ALL CATALOG TABLES TESTS PASSED")


if __name__ == "__main__":
    main()
