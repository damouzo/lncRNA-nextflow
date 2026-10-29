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
        "--extreme-basemean", "1000", "--genome-build", "GRCh38",
        "--out", os.path.join(out_dir, "de_full.tsv"),
        "--significant-out", os.path.join(out_dir, "de_significant.tsv"),
        "--lncrna-out", os.path.join(out_dir, "de_lncrna_full.tsv"),
        "--summary-out", os.path.join(out_dir, "summary.tsv"),
    ], capture_output=True, text=True)


def main():
    d = tempfile.mkdtemp()
    write(os.path.join(d, "catalog.tsv"),
          "gene_id\torigin\tclass_code\ttier_static\ttier_final\tsingle_exon_low_confidence\tantisense_status\t"
          "gene_biotype\tcatalog_id\tchromosome\tstart\tend\tstrand\tlocus\n"
          "G1\tnovel\tu\tTRUE\tTRUE\tFALSE\tnot_evaluated\tNA\tLNC000001\tKI270733.1\t100\t200\t+\tKI270733.1:100-200(+)\n"
          "G2\tnovel\tu\tFALSE\tFALSE\tFALSE\tnot_evaluated\tNA\tLNC000002\t1\t500\t600\t-\t1:500-600(-)\n"
          "G3\tannotated_lncrna\tNA\tNA\tFALSE\tNA\tnot_evaluated\tlncRNA\tNA\t4\t1000\t2000\t+\t4:1000-2000(+)\n")

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
            check(all(c in full[0] for c in ("stat", "lfc_unshrunken", "shrink_type",
                                             "gene_name", "gene_biotype", "catalog_id",
                                             "lnc_position_class")),
                  "%s__%s: unshrunken stats + extra annotation columns present"
                  % (group, ct))
            # 4a.3: both significance definitions, plus the strict alias.
            check(by["G1"]["significant_padj"] == "TRUE" and by["G1"]["significant_padj_lfc1"] == "TRUE",
                  "%s__%s: G1 significant under both definitions" % (group, ct))
            check(by["G2"]["significant_padj"] == "TRUE" and by["G2"]["significant_padj_lfc1"] == "FALSE",
                  "%s__%s: G2 significant by padj only (LFC not met)" % (group, ct))
            check(by["G1"]["significant"] == by["G1"]["significant_padj_lfc1"],
                  "%s__%s: `significant` == significant_padj_lfc1" % (group, ct))
            # 4b.1/4b.4: positions and gene_group are joined from the catalog.
            check(by["G3"]["chr"] == "4" and by["G3"]["gene_group"] == "annotated_lncRNA",
                  "%s__%s: annotated lncRNA carries position + gene_group" % (group, ct))
            check(by["G1"]["gene_group"] == "novel_lncRNA" and by["G1"]["chr"] == "KI270733.1",
                  "%s__%s: novel gene_group + contig position" % (group, ct))
            # 4b.1 (decisión b): one stable id column.
            check(by["G1"]["gene_id_stable"] == "LNC000001" and by["G3"]["gene_id_stable"] == "G3",
                  "%s__%s: gene_id_stable = catalog_id for novel, gene_id otherwise" % (group, ct))
            check(by["G1"]["locus"] == "KI270733.1:100-200(+)",
                  "%s__%s: locus joined from the catalog" % (group, ct))
            check(all(r["genome_build"] == "GRCh38" for r in full),
                  "%s__%s: genome_build recorded on every row" % (group, ct))
            # 4a.4: non-canonical contig flag + extreme baseMean flag.
            check(by["G1"]["contig_type"] == "noncanonical" and by["G2"]["contig_type"] == "canonical",
                  "%s__%s: contig_type flags KI270733.1" % (group, ct))
            check(by["G2"]["baseMean_extreme"] == "FALSE",
                  "%s__%s: baseMean below threshold is not flagged" % (group, ct))
            # 4b.4: de_lncrna_full keeps only lncRNAs, unfiltered by significance.
            lnc = read_tsv(os.path.join(out_dir, "de_lncrna_full.tsv"))
            check({r["gene_id"] for r in lnc} == {"G1", "G2", "G3"},
                  "%s__%s: de_lncrna_full holds all lncRNAs (incl. non-significant)"
                  % (group, ct))
            check(all(r["gene_group"] in ("annotated_lncRNA", "novel_lncRNA") for r in lnc),
                  "%s__%s: de_lncrna_full has no protein_coding/other rows" % (group, ct))

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
