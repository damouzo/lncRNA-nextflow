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

    # 1.5: categorical tier {1a,1b,2} for novel genes, NA for reference rows.
    g1 = [g for g in genes if g["gene_id"] == "G1"][0]
    check(g1["tier"] in {"1a", "1b", "2"} and g1["tier"] == "1a",
          "novel gene carries the categorical tier from novel_type")
    check(g2["tier"] == "NA", "reference gene tier is NA (not the literal None)")
    check(all(g["tier"] in {"1a", "1b", "2", "NA"} for g in genes),
          "tier only takes 1a/1b/2/NA")

    _, novel = B.read_tsv(os.path.join(out, "novel_lncrnas.tsv"))
    check([r["gene_id"] for r in novel] == ["G1"],
          "novel_lncrnas.tsv is the origin==novel view of catalog_genes")

    n_tx = B.build_catalog_tx(os.path.join(d, "tier_final_transcripts.tsv"), out)
    check(n_tx == 1, "catalog_tx: 1 transcript")
    _, tx = B.read_tsv(os.path.join(out, "catalog_tx.tsv"))
    check(all(c in tx[0] for c in ("catalog_id", "catalog_transcript_id",
                                   "mstrg_transcript_id", "closest_reference_transcript")),
          "catalog_tx carries the stable-id columns")
    check(tx[0]["catalog_id"] == "NA",
          "catalog_tx stable ids are NA without a catalog GTF")

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

    # Freeze class gate: a k transcript that passes coding + recurrence is not frozen.
    write(os.path.join(d, "gffcmp_k.gtf"),
          'chr1\ts\ttranscript\t1\t300\t.\t+\t.\ttranscript_id "T1"; class_code "u";\n'
          'chr1\ts\texon\t1\t300\t.\t+\t.\ttranscript_id "T1";\n'
          'chr1\ts\ttranscript\t1\t400\t.\t+\t.\ttranscript_id "T5"; class_code "k";\n'
          'chr1\ts\texon\t1\t400\t.\t+\t.\ttranscript_id "T5";\n'
          'chr1\ts\ttranscript\t1\t400\t.\t+\t.\ttranscript_id "T6"; class_code "o";\n'
          'chr1\ts\texon\t1\t400\t.\t+\t.\ttranscript_id "T6";\n')
    write(os.path.join(d, "recurrence_k.tsv"),
          "ID\tconsensus_noncoding\tpasses_recurrence\n"
          "T1\tTRUE\tTRUE\nT5\tTRUE\tTRUE\nT6\tTRUE\tTRUE\n")
    _, rows_k = B.build_candidates_all(os.path.join(d, "gffcmp_k.gtf"),
                                       os.path.join(d, "recurrence_k.tsv"),
                                       os.path.join(d, "frozen.gtf"), 200, out,
                                       novel_classes={"u", "x", "i", "j"})
    bk = {r["transcript_id"]: r for r in rows_k}
    check(bk["T1"]["filter_status"] == "frozen", "novel-class u still frozen")
    check(bk["T5"]["filter_status"] == "filtered" and bk["T5"]["filter_reason"] == "not_novel_class"
          and bk["T5"]["frozen"] == "FALSE",
          "k passing every filter is not_novel_class, never frozen")
    check(bk["T6"]["filter_reason"] == "not_novel_class", "o is not_novel_class")
    _, rows_nogate = B.build_candidates_all(os.path.join(d, "gffcmp_k.gtf"),
                                            os.path.join(d, "recurrence_k.tsv"),
                                            os.path.join(d, "frozen.gtf"), 200, out)
    check({r["transcript_id"]: r for r in rows_nogate}["T5"]["filter_reason"] == "",
          "without novel_classes the class gate is off (back-compat)")

    try:
        B.build_candidates_all(os.path.join(d, "gffcmp.annotated.gtf"),
                               os.path.join(d, "recurrence_filtered.tsv"),
                               os.path.join(d, "frozen.gtf"), 200, out,
                               os.path.join(d, "length_filtered_short.gtf"))
        check(False, "parity mismatch must raise SystemExit")
    except SystemExit:
        check(True, "parity mismatch raises SystemExit")


def test_catalog_tx_ids():
    # catalog_tx.tsv resolves the LNC ids + cmp_ref from the ASSIGN_CATALOG_IDS GTF.
    d = tempfile.mkdtemp()
    write(os.path.join(d, "tier_final_transcripts.tsv"),
          "transcript_id\tgene_id\tclass_code\tnovel_type\ttier_static\ttier_final\n"
          "MSTRG.351__ENST00000584329\tMSTRG.351\tk\t2\tFALSE\tFALSE\n")
    write(os.path.join(d, "catalog.gtf"),
          'chr1\ts\ttranscript\t1\t500\t.\t+\t.\ttranscript_id "LNC000007.1"; gene_id "LNC000007"; '
          'cmp_ref "ENST00000584329"; class_code "k"; mstrg_gene_id "MSTRG.351"; '
          'mstrg_transcript_id "MSTRG.351__ENST00000584329"; catalog_id "LNC000007";\n')
    out = os.path.join(d, "out")
    cat_map = B.parse_catalog_gtf_ids(os.path.join(d, "catalog.gtf"))
    B.build_catalog_tx(os.path.join(d, "tier_final_transcripts.tsv"), out, cat_map)
    _, tx = B.read_tsv(os.path.join(out, "catalog_tx.tsv"))
    r = tx[0]
    check(r["catalog_id"] == "LNC000007", "catalog_tx resolves gene catalog_id")
    check(r["catalog_transcript_id"] == "LNC000007.1", "catalog_tx resolves transcript catalog id")
    check(r["mstrg_transcript_id"] == "MSTRG.351__ENST00000584329", "catalog_tx keeps the mstrg id")
    check(r["closest_reference_transcript"] == "ENST00000584329",
          "catalog_tx carries the closest reference transcript (cmp_ref)")


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
    # 4b.2b: gene_group + a built locus for every row.
    check(by["G1"]["gene_group"] == "novel_lncRNA", "novel row gene_group")
    check(by["CG1"]["gene_group"] == "protein_coding", "coding row gene_group")
    check(by["CG1"]["locus"] == "1:1000-2000(+)", "reference row gets a locus built from coords")
    check(by["G1"]["locus"] == "1:1500-1600(-)", "novel row locus built from coords")
    check(all(g["gene_group"] in ("novel_lncRNA", "annotated_lncRNA", "protein_coding", "other")
              for g in genes), "gene_group only takes the four allowed values")


def test_tier_consistency_abort():
    # H26 / 1.7: BUILD_GENE_CATALOG (R) vs COMPUTE_TIER_FINAL (Python) drift -> abort.
    d = tempfile.mkdtemp()
    write(os.path.join(d, "gene_catalog.tsv"),
          "gene_id\torigin\tclass_code\tn_tier1_transcripts\ttier_static\tnovel_type\n"
          "G1\tnovel\tu\t2\tTRUE\t1a\n")
    write(os.path.join(d, "tier_final_genes.tsv"),
          "gene_id\ttier_static\tn_tier1_transcripts\tantisense_status\ttier_final\n"
          "G1\tFALSE\t2\tevaluated\tFALSE\n")
    out = os.path.join(d, "out")
    tier_final = B.index_unique(os.path.join(d, "tier_final_genes.tsv"), "gene_id")
    try:
        B.build_catalog_genes(os.path.join(d, "gene_catalog.tsv"), tier_final, out)
        check(False, "tier_static mismatch must raise SystemExit")
    except SystemExit:
        check(True, "tier_static mismatch raises SystemExit (H26 abort)")


def test_gene_annotation_joined():
    # catalog_tx.tsv and candidates_all.tsv must carry the same gene-level
    # classification columns as catalog_genes.tsv.
    d = tempfile.mkdtemp()
    make_fixture(d)
    out = os.path.join(d, "out")
    tier_final = B.index_unique(os.path.join(d, "tier_final_genes.tsv"), "gene_id")
    B.build_catalog_genes(os.path.join(d, "gene_catalog.tsv"), tier_final, out)
    gene_ann = B.load_gene_annotation(os.path.join(out, "catalog_genes.tsv"))
    B.build_catalog_tx(os.path.join(d, "tier_final_transcripts.tsv"), out, None, gene_ann)
    _, tx = B.read_tsv(os.path.join(out, "catalog_tx.tsv"))
    check(tx[0]["origin"] == "novel", "catalog_tx carries the gene origin")
    check(tx[0]["gene_group"] == "novel_lncRNA", "catalog_tx carries gene_group")
    check("lnc_position_class" in tx[0], "catalog_tx carries lnc_position_class")

    _, rows = B.build_candidates_all(os.path.join(d, "gffcmp.annotated.gtf"),
                                     os.path.join(d, "recurrence_filtered.tsv"),
                                     os.path.join(d, "frozen.gtf"), 200, out,
                                     None, gene_ann)
    by = {r["transcript_id"]: r for r in rows}
    check(all(k in by["T1"] for k in ("gene_id", "catalog_id", "origin",
                                      "gene_group", "tier", "lnc_position_class")),
          "candidates_all carries gene_id + the annotation columns")


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
    for f in ("catalog_genes.tsv", "catalog_tx.tsv", "candidates_all.tsv", "novel_lncrnas.tsv"):
        check(os.path.exists(os.path.join(out, f)), "CLI wrote %s" % f)


def main():
    for fn in (test_build, test_catalog_tx_ids, test_position_class_wired,
               test_gene_annotation_joined, test_tier_consistency_abort, test_cli):
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
