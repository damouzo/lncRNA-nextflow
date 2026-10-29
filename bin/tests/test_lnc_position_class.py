#!/usr/bin/env python3
"""F5.4 (lnc_position_class) tests. Stdlib only.

Run: python3 bin/tests/test_lnc_position_class.py
"""

import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.dirname(HERE)
sys.path.insert(0, BIN)

import lnc_position_class as L  # noqa: E402

FAILS = []


def check(cond, msg):
    print(("  ok  " if cond else "  FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def coding(gid, chrom, start, end, strand, biotype="protein_coding",
           origin="reference_other"):
    return {"gene_id": gid, "origin": origin, "gene_biotype": biotype,
            "chromosome": chrom, "start": str(start), "end": str(end), "strand": strand}


def lnc(gid, chrom, start, end, strand, origin="novel"):
    return {"gene_id": gid, "origin": origin, "gene_biotype": "NA",
            "chromosome": chrom, "start": str(start), "end": str(end), "strand": strand}


def classify(rows, gid):
    L.annotate_rows(rows)
    return next(r for r in rows if r["gene_id"] == gid)


def test_basic_classes():
    rows = [
        coding("CG1", "chr1", 1000, 2000, "+"),
        lnc("A", "chr1", 1500, 1600, "-"),   # opposite strand -> antisense
        lnc("I", "chr1", 1200, 1300, "+"),   # contained, same strand -> intronic
        lnc("O", "chr1", 1900, 2500, "+"),   # partial, same strand -> overlapping
        lnc("G", "chr1", 9000, 9100, "+"),   # no overlap -> no_coding_overlap
    ]
    check(classify(rows, "A")["lnc_position_class"] == "antisense", "opposite strand -> antisense")
    check(classify(rows, "I")["lnc_position_class"] == "intronic", "contained same strand -> intronic")
    check(classify(rows, "O")["lnc_position_class"] == "overlapping", "partial same strand -> overlapping")
    g = classify(rows, "G")
    check(g["lnc_position_class"] == "no_coding_overlap", "no overlap -> no_coding_overlap")
    check(g["lnc_position_gene_id"] == "NA" and g["lnc_position_overlap_bp"] == "0"
          and g["lnc_position_n_genes"] == "0" and g["lnc_position_genes"] == "",
          "no_coding_overlap carries NA/0/empty audit fields")


def test_priority_and_tiebreak():
    # antisense beats a larger intronic overlap
    rows = [
        coding("CG_int", "chr1", 50, 300, "+"),
        coding("CG_anti", "chr1", 150, 250, "-"),
        lnc("A", "chr1", 100, 200, "+"),
    ]
    check(classify(rows, "A")["lnc_position_class"] == "antisense",
          "antisense priority beats larger intronic overlap")

    # same class: larger overlap_bp wins
    rows = [
        coding("CG_A", "chr1", 100, 200, "+"),
        coding("CG_B", "chr1", 150, 400, "+"),
        lnc("O", "chr1", 100, 300, "+"),
    ]
    r = classify(rows, "O")
    check(r["lnc_position_gene_id"] == "CG_B", "same class -> larger overlap_bp wins")
    check(r["lnc_position_overlap_bp"] == "151", "overlap_bp recorded for the winner")

    # identical overlap: smaller gene_id wins
    rows = [
        coding("CG_Y", "chr1", 100, 200, "+"),
        coding("CG_X", "chr1", 100, 200, "+"),
        lnc("T", "chr1", 100, 150, "+"),
    ]
    check(classify(rows, "T")["lnc_position_gene_id"] == "CG_X", "tie -> smaller gene_id wins")


def test_multiple_overlaps_audited():
    rows = [
        coding("CG1", "chr1", 100, 500, "+"),
        coding("CG2", "chr2", 100, 500, "+"),
        coding("CG3", "chr1", 600, 900, "+"),
        lnc("L", "chr1", 400, 700, "-"),
    ]
    r = classify(rows, "L")
    check(r["lnc_position_n_genes"] == "2", "n_genes counts all overlapping coding genes")
    check(r["lnc_position_genes"] == "CG1,CG3", "overlapping gene list is sorted")
    check(r["lnc_position_class"] == "antisense", "opposite strand wins over the sense partner")


def test_lower_bound_skips_distant_genes():
    rows = [
        coding("CG_far", "chr1", 100, 200, "+"),
        coding("CG_near", "chr1", 3000, 4000, "+"),
        lnc("L", "chr1", 3500, 3600, "+"),
    ]
    r = classify(rows, "L")
    check(r["lnc_position_n_genes"] == "1" and r["lnc_position_gene_id"] == "CG_near",
          "lower bound excludes distant non-overlapping genes")


def test_origin_filtering():
    rows = [
        coding("CG1", "chr1", 1000, 2000, "+"),
        coding("CGalt", "chr1", 1000, 2000, "+", biotype="lncRNA"),   # not coding -> ignored
        lnc("N", "chr1", 1500, 1600, "+"),                             # intronic to CG1
        lnc("AL", "chr1", 1500, 1600, "-", origin="annotated_lncrna"),  # annotated lncRNA classified too
        coding("R", "chr1", 3000, 4000, "-"),                          # coding row itself -> NA
    ]
    L.annotate_rows(rows)
    by = {r["gene_id"]: r for r in rows}
    check(by["N"]["lnc_position_class"] == "intronic", "novel lncRNA classified")
    check(by["AL"]["lnc_position_class"] == "antisense", "annotated_lncrna classified")
    check(by["AL"]["lnc_position_gene_id"] == "CG1",
          "non-coding reference_other is not used as a coding partner")
    check(by["R"]["lnc_position_class"] == "NA", "reference_other coding row gets NA")
    check(by["CGalt"]["lnc_position_class"] == "NA",
          "reference_other non-coding row gets NA")


def test_missing_coordinates():
    rows = [
        coding("CG1", "chr1", 1000, 2000, "+"),
        lnc("N", "chr1", "NA", "NA", "+"),
    ]
    check(classify(rows, "N")["lnc_position_class"] == "NA",
          "lncRNA without coordinates gets NA")


def test_cli():
    with tempfile.TemporaryDirectory() as d:
        cat = os.path.join(d, "catalog_genes.tsv")
        with open(cat, "w") as fh:
            fh.write("gene_id\torigin\tgene_biotype\tchromosome\tstart\tend\tstrand\n")
            fh.write("CG1\treference_other\tprotein_coding\tchr1\t1000\t2000\t+\n")
            fh.write("N\tnovel\tNA\tchr1\t1500\t1600\t-\n")
        out = os.path.join(d, "out.tsv")
        rc = subprocess.run([sys.executable, os.path.join(BIN, "lnc_position_class.py"),
                             "--catalog", cat, "--out", out]).returncode
        check(rc == 0, "CLI exits 0")
        with open(out) as fh:
            lines = fh.read().splitlines()
        check("lnc_position_class" in lines[0], "CLI appends the position columns")
        check(any("antisense" in l for l in lines[1:]), "CLI classifies the lncRNA")


def test_overlap_with_non_coding_annotated_is_not_intergenic():
    # H13 / 1.6: only protein-coding genes are compared, so an overlap with an
    # annotated lncRNA or a pseudogene must not be reported as `intergenic`.
    rows = [
        lnc("L", "chr1", 1000, 2000, "+", origin="novel"),
        lnc("ANN", "chr1", 1200, 1300, "+", origin="annotated_lncrna"),
        coding("PSEUDO", "chr1", 1400, 1500, "+", biotype="pseudogene",
               origin="reference_other"),
    ]
    r = classify(rows, "L")
    check(r["lnc_position_class"] == "no_coding_overlap",
          "overlap with lncRNA/pseudogene is no_coding_overlap, not intergenic")


def main():
    for t in (test_basic_classes, test_priority_and_tiebreak, test_multiple_overlaps_audited,
              test_lower_bound_skips_distant_genes, test_origin_filtering,
              test_missing_coordinates,
              test_overlap_with_non_coding_annotated_is_not_intergenic, test_cli):
        t()
    print()
    if FAILS:
        print("LNCRNA POSITION TESTS FAILED: %d" % len(FAILS))
        return 1
    print("ALL LNCRNA POSITION TESTS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
