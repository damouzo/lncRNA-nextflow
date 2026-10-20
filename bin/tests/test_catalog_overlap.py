#!/usr/bin/env python3
"""F2.1 reciprocal-overlap tests on synthetic intervals. Stdlib only.

Run: python3 bin/tests/test_catalog_overlap.py
"""

import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import catalog_overlap as O  # noqa: E402

FAILS = []


def check(cond, msg):
    print(("  ok  " if cond else "  FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def test_reciprocal():
    print("reciprocal")
    # small candidate inside a large feature -> candidate side is 1.0, feature side tiny
    st, perm, ov, fc, fe = O.reciprocal([(120, 130)], [(100, 200)])
    check(perm and not st and fe < 0.5, "small candidate inside large feature is not strict")
    # reciprocal both >= 0.5
    st, perm, ov, fc, fe = O.reciprocal([(100, 200)], [(120, 180)])
    check(st and fc >= 0.5 and fe >= 0.5, "mutual >=0.5 is strict")
    # opposite strand handled by caller; intervals alone still overlap
    st, perm, _, _, _ = O.reciprocal([(100, 200)], [(300, 400)])
    check(not perm, "no overlap -> not permissive")


def _cand(exons, strand="+"):
    return {"tid": {"chrom": "1", "strand": strand, "exons": O.C.merge_intervals(exons)}}


def test_not_summed():
    print("features are not summed")
    cand = _cand([(100, 200)])
    ext = [{"chrom": "1", "strand": "+", "name": "A", "exons": [(100, 140)], "span_only": False},
           {"chrom": "1", "strand": "+", "name": "B", "exons": [(161, 200)], "span_only": False}]
    r = O.annotate_source(cand, ext, "x", 0.5)["tid"]
    check(r["strict"] is False, "two half-features that jointly cover are not a strict hit")
    check(r["permissive"] is True, "but it is permissive")


def test_order_and_best():
    print("order irrelevant + best pair")
    cand = _cand([(100, 200)])
    ext = [{"chrom": "1", "strand": "+", "name": "A", "exons": [(100, 140)], "span_only": False},
           {"chrom": "1", "strand": "+", "name": "B", "exons": [(130, 190)], "span_only": False}]
    r1 = O.annotate_source(cand, ext, "x", 0.5)["tid"]
    r2 = O.annotate_source(cand, list(reversed(ext)), "x", 0.5)["tid"]
    check(r1["ext_id"] == r2["ext_id"], "order of features does not change the best pair")
    check(r1["ext_id"] == "B", "best pair chosen by max min(frac_candidate, frac_external)")


def test_strand():
    print("strand")
    cand = _cand([(100, 200)], strand="+")
    ext = [{"chrom": "1", "strand": "-", "name": "A", "exons": [(100, 200)], "span_only": False}]
    r = O.annotate_source(cand, ext, "x", 0.5)["tid"]
    check(r["permissive"] is False and r["strict"] is False, "opposite strand -> no hit")


def test_parse_bed():
    print("parse BED6 / BED12")
    d = tempfile.mkdtemp()
    bed6 = os.path.join(d, "a.bed")
    with open(bed6, "w") as fh:
        fh.write("1\t99\t200\tfeatA\t0\t+\n")
    fa = O.parse_bed(bed6)[0]
    check(fa["span_only"] is True and fa["exons"] == [(100, 200)], "BED6 -> span, marked span_only")
    bed12 = os.path.join(d, "b.bed")
    with open(bed12, "w") as fh:
        fh.write("1\t99\t200\tfeatB\t0\t+\t99\t200\t0\t2\t50,50,\t0,51,\n")
    fb = O.parse_bed(bed12)[0]
    check(fb["span_only"] is False and fb["exons"] == [(100, 149), (151, 200)], "BED12 -> exon blocks")


def test_gtf_parser():
    print("parse external GTF exons")
    d = tempfile.mkdtemp()
    gtf = os.path.join(d, "ext.gtf")
    with open(gtf, "w") as fh:
        fh.write('1\tX\ttranscript\t100\t200\t.\t+\t.\tgene_id "g"; transcript_id "t";\n')
        fh.write('1\tX\texon\t100\t149\t.\t+\t.\tgene_id "g"; transcript_id "t";\n')
        fh.write('1\tX\texon\t151\t200\t.\t+\t.\tgene_id "g"; transcript_id "t";\n')
    fx = O.parse_gtf_exons(gtf)
    check(len(fx) == 1 and fx[0]["exons"] == [(100, 149), (151, 200)], "GTF exons parsed per transcript")


def test_consensus_annotation():
    print("consensus annotation (drop-in for CATALOG_OVERLAP)")
    d = tempfile.mkdtemp()
    cons = os.path.join(d, "consensus.tsv")
    with open(cons, "w") as fh:
        fh.write("ID\tclass_code\nT1\ti\nT2\tx\n")
    gtf = os.path.join(d, "consensus.gtf")
    with open(gtf, "w") as fh:
        fh.write('1\tS\ttranscript\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "T1";\n')
        fh.write('1\tS\texon\t100\t120\t.\t+\t.\tgene_id "g1"; transcript_id "T1";\n')
        fh.write('1\tS\texon\t180\t200\t.\t+\t.\tgene_id "g1"; transcript_id "T1";\n')
        fh.write('1\tS\ttranscript\t100\t200\t.\t+\t.\tgene_id "g2"; transcript_id "T2";\n')
        fh.write('1\tS\texon\t100\t200\t.\t+\t.\tgene_id "g2"; transcript_id "T2";\n')
    bed = os.path.join(d, "lncipedia.bed")
    with open(bed, "w") as fh:
        fh.write("1\t99\t200\tLNCT2\t0\t+\n")
    header, rows = O.C.read_tsv(cons)
    cands = O.exons_by_tx(gtf)
    nh, out = O.annotate_consensus(header, rows, cands, {"lncipedia": bed}, 0.5)
    by = {r["ID"]: r for r in out}
    check("known_catalog_hit" in nh and "lncipedia_pct" in nh, "columns appended")
    check(by["T2"]["known_catalog_hit"] == "TRUE" and by["T2"]["lncipedia_id"] == "LNCT2",
          "full-cover transcript is a strict hit")
    check(by["T1"]["known_catalog_hit"] == "FALSE",
          "candidate whose exon union is only partly covered is not strict")
    check(by["T2"]["catalog_sources"] == "lncipedia", "catalog_sources lists the source")


def main():
    for fn in (test_reciprocal, test_not_summed, test_order_and_best, test_strand,
               test_parse_bed, test_gtf_parser, test_consensus_annotation):
        fn()
    print()
    if FAILS:
        print("FAILED: %d" % len(FAILS))
        sys.exit(1)
    print("ALL CATALOG OVERLAP TESTS PASSED")


if __name__ == "__main__":
    main()
