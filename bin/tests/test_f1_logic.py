#!/usr/bin/env python3
"""Pure-logic tests for the F1 helpers. No R, no samtools, no Nextflow.

Run: python3 bin/tests/test_f1_logic.py
Exits non-zero on the first failure. Uses only the standard library.
"""

import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import f1_common as C          # noqa: E402
import f1_regions as R         # noqa: E402
import f1_metrics as M         # noqa: E402
import f1_thresholds as T      # noqa: E402
import f1_tier_final as F      # noqa: E402

FAILS = []


def check(cond, msg):
    if cond:
        print("  ok  %s" % msg)
    else:
        print("  FAIL %s" % msg)
        FAILS.append(msg)


def test_intervals():
    print("intervals")
    check(C.merge_intervals([(1, 5), (6, 10), (20, 25)]) == [(1, 10), (20, 25)], "merge adjacent")
    check(C.merge_intervals([(1, 10), (3, 5)]) == [(1, 10)], "merge contained")
    check(C.overlaps((1, 5), (5, 9)), "overlap touching at boundary")
    check(not C.overlaps((1, 4), (5, 9)), "no overlap")
    check(C.subtract_intervals((10, 20), [(12, 14)]) == [(10, 11), (15, 20)], "subtract middle")
    check(C.subtract_intervals((10, 20), [(1, 12), (18, 30)]) == [(13, 17)], "subtract clipped")
    check(C.subtract_intervals((10, 20), [(10, 20)]) == [], "subtract all")
    check(C.subtract_intervals((10, 20), [(11, 12), (12, 13)]) == [(10, 10), (14, 20)], "subtract adjacent masks")
    check(C.interval_intersection((1, 10), (5, 15)) == (5, 10), "intersection")
    check(C.interval_intersection((1, 4), (5, 9)) is None, "intersection none")


def test_stats():
    print("stats")
    check(C.median([3, 1, 2]) == 2.0, "median odd")
    check(C.median([1, 2, 3, 4]) == 2.5, "median even")
    check(C.median([]) is None, "median empty")
    check(abs(C.percentile([0, 10], 50) - 5.0) < 1e-9, "percentile 50")
    check(C.percentile([5], 10) == 5.0, "percentile single")
    check(T.rejection_rate([1, 2, 3, 4], 3) == 0.5, "rejection rate")
    check(T.rejection_rate([], 3) is None, "rejection rate empty")
    check(abs(T.flag_rate([0.1, 0.5, 0.4], 0.32) - 1.0 / 3) < 1e-9, "flag rate")
    check(T.flag_rate([], 0.3) is None, "flag rate empty")
    rho = C.spearman([1, 2, 3, 4, 5], [2, 4, 6, 8, 10])
    check(rho is not None and abs(rho - 1.0) < 1e-9, "spearman perfect")
    check(C.spearman([1, 2], [2, 1]) is None, "spearman n<3")


def test_attrs():
    print("gtf attributes")
    a = C.parse_attributes('gene_id "G1"; transcript_id "T1"; class_code "x";')
    check(a["gene_id"] == "G1" and a["class_code"] == "x", "parse attributes")
    check(C.strip_version("ENST0001.5") == "ENST0001", "strip version")
    check(C.ref_gene_of({"ref_gene_id": "G1.2"}) == "G1", "ref_gene_of prefers ref_gene_id")
    check(C.ref_gene_of({"cmp_ref_gene": "G2"}) == "G2", "ref_gene_of falls back to cmp_ref_gene")
    check(C.ref_gene_of({}) is None, "ref_gene_of empty -> None")
    tx2gene = {"T1": "G1"}
    check(C.host_gene_of({"ref_gene_id": "G9.3"}, tx2gene) == "G9", "host_gene_of prefers ref_gene_id")
    check(C.host_gene_of({"cmp_ref": "T1.5"}, tx2gene) == "G1", "host_gene_of maps cmp_ref via tx2gene")
    check(C.host_gene_of({"cmp_ref_gene": "SYM"}, tx2gene) is None,
          "host_gene_of ignores cmp_ref_gene (symbol)")


def test_regions_pure():
    print("regions pure")
    check(R.genomic_flanks(100, 200, 50, 1000) == [("left", (50, 99)), ("right", (201, 250))],
          "genomic flanks")
    check(R.genomic_flanks(10, 20, 50, 1000) == [("left", (1, 9)), ("right", (21, 70))],
          "flank clipped at start")
    check(R.mask_pieces((100, 200), [(120, 130)], 50) == [(100, 119), (131, 200)], "mask keep")
    check(R.mask_pieces((100, 200), [(100, 200)], 50) is None, "mask all -> None")
    check(R.mask_pieces((100, 110), [], 50) is None, "too short -> None")
    check(R.host_intron(150, 160, [(100, 120), (200, 300)]) == (121, 199), "host intron")
    check(R.host_intron(10, 20, [(100, 120)]) is None, "host intron none")
    check(R.antisense_region([(100, 200)], [(50, 150)]) == [(100, 150)], "antisense region")
    check(R.antisense_region([(100, 110)], [(200, 210)]) is None, "antisense none")
    t1a, t1, pb = {"u", "x", "i"}, {"u", "x", "i", "j", "o"}, {"lncRNA"}
    check(C.tier_of("u", None, t1a, t1, pb) == "1a", "tier 1a")
    check(C.tier_of("j", "lncRNA", t1a, t1, pb) == "1b", "tier 1b lncRNA parent")
    check(C.tier_of("j", "protein_coding", t1a, t1, pb) is None, "j coding parent -> tier2")
    check(C.tier_of("o", "processed_pseudogene", t1a, t1, pb) is None, "o pseudogene parent -> tier2")
    check(C.tier_of("k", "lncRNA", t1a, t1, pb) is None, "k not tier1")
    check(C.tier_of("j", None, t1a, t1, pb) is None, "j unknown parent -> tier2")


def test_regions_integration():
    print("regions integration")
    d = tempfile.mkdtemp()
    ref = os.path.join(d, "ref.gtf")
    frozen = os.path.join(d, "frozen.gtf")
    rec = os.path.join(d, "rec.tsv")
    with open(ref, "w") as fh:
        fh.write('1\tX\tgene\t1000\t1500\t.\t+\t.\tgene_id "G1"; gene_biotype "protein_coding";\n')
        fh.write('1\tX\ttranscript\t1000\t1500\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n')
        fh.write('1\tX\texon\t1000\t1100\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n')
        fh.write('1\tX\texon\t1400\t1500\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n')
        fh.write('1\tX\tgene\t8000\t8200\t.\t-\t.\tgene_id "G2"; gene_biotype "lncRNA";\n')
        fh.write('1\tX\ttranscript\t8000\t8200\t.\t-\t.\tgene_id "G2"; transcript_id "T2";\n')
        fh.write('1\tX\texon\t8000\t8200\t.\t-\t.\tgene_id "G2"; transcript_id "T2";\n')
    with open(frozen, "w") as fh:
        fh.write('1\tS\ttranscript\t3000\t3200\t.\t+\t.\tgene_id "M1"; transcript_id "M1.1"; class_code "u"; ref_gene_id "G1";\n')
        fh.write('1\tS\texon\t3000\t3200\t.\t+\t.\tgene_id "M1"; transcript_id "M1.1";\n')
        fh.write('1\tS\ttranscript\t1000\t1100\t.\t-\t.\tgene_id "M2"; transcript_id "M2.1"; class_code "x"; ref_gene_id "G1";\n')
        fh.write('1\tS\texon\t1000\t1100\t.\t-\t.\tgene_id "M2"; transcript_id "M2.1";\n')
    with open(rec, "w") as fh:
        fh.write("ID\tn_reps_expressed\tpasses_recurrence\tclass_code\n")
        fh.write("T1\t5\tTRUE\t=\n")
        fh.write("T2\t5\tTRUE\t=\n")
    with open(os.path.join(d, "done"), "w") as fh:
        fh.write("")
    os.makedirs(os.path.join(d, "out"), exist_ok=True)

    class A:
        pass
    a = A()
    a.tier1_classes = ["u", "x", "i", "j", "o"]
    a.tier1a_classes = ["u", "x", "i"]
    a.tier1_parent_biotypes = ["lncRNA"]
    a.annotated_lncrna_biotypes = ["lncRNA"]
    a.flank_bp = 500
    a.min_flank_len = 50
    a.cap = 100
    a.seed = 1
    a.expr_min_reps = 3
    reference = C.load_reference(ref, ["lncRNA"])
    frozen_g = C.load_gtf(frozen, {"transcript", "exon"})
    R.build(reference, frozen_g, rec, os.path.join(d, "out"), a)
    _, regions = C.read_tsv(os.path.join(d, "out", "regions.tsv"))
    ids = {r["element_id"] for r in regions}
    check("cand:M1.1" in ids, "candidate u present")
    check("cand_antisense:M2.1" in ids, "candidate x antisense present")
    check(any(r["set"] == "pos_flank" for r in regions), "pos_flank control present")
    check(any(r["set"] == "neg_flank" for r in regions), "neg_flank control present")
    check(any(r["set"] == "null_antisense" for r in regions), "null_antisense control present")
    # determinism: rebuilding yields identical regions.tsv
    R.build(reference, frozen_g, rec, os.path.join(d, "out2"), a)
    with open(os.path.join(d, "out", "regions.tsv")) as f1, \
         open(os.path.join(d, "out2", "regions.tsv")) as f2:
        check(f1.read() == f2.read(), "regions deterministic")


def _fc(path, rows):
    with open(path, "w") as fh:
        fh.write("# Program:featureCounts\n")
        fh.write("Geneid\tChr\tStart\tEnd\tStrand\tLength\tbamA.bam\tbamB.bam\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


def test_metrics():
    print("metrics")
    d = tempfile.mkdtemp()
    s2 = os.path.join(d, "s2.tsv")
    s1 = os.path.join(d, "s1.tsv")
    # element e1: candidate 100 reads in A, 200 in B; flanks 10 each
    rows2 = [
        ("e1|candidate", "1", 100, 199, "+", 100, 100, 200),
        ("e1|left", "1", 1, 50, "+", 50, 10, 10),
        ("e1|right", "1", 250, 299, "+", 50, 10, 10),
        ("a1|region", "1", 100, 199, "+", 100, 50, 100),
    ]
    rows1 = [
        ("e1|candidate", "1", 100, 199, "+", 100, 10, 20),
        ("a1|region", "1", 100, 199, "+", 100, 5, 10),
    ]
    _fc(s2, rows2)
    _fc(s1, rows1)
    regions = os.path.join(d, "regions.tsv")
    C.write_tsv([("e1", "flank", "candidate", "1", "+", ""),
                 ("a1", "antisense", "candidate", "1", "+", "")],
                regions, ["element_id", "kind", "set", "chrom", "strand", "note"])
    coldata = os.path.join(d, "coldata.tsv")
    with open(coldata, "w") as fh:
        fh.write("S1\tgrp1\tbamA.bam\nS2\tgrp1\tbamB.bam\n")

    _, reg = C.read_tsv(regions)
    regmap = {r["element_id"]: r for r in reg}
    bams = M.read_coldata(coldata)
    meta, c2, cols = M.parse_featurecounts(s2)
    _, c1, _ = M.parse_featurecounts(s1)
    fm = M.flank_metrics(meta, c2, c1, cols, bams, regmap, min_reads=1)
    # S1: cand density 100/100=1 ; flank density 10/50=0.2 -> ratio 5
    # S2: cand density 2 ; flank density 0.2 -> ratio 10 ; opp ratio S1=0.1
    row = [r for r in fm if r[0] == "e1"][0]
    check(row[0] == "e1" and row[2] == "grp1", "flank row keyed")
    check(abs(float(row[3]) - 7.5) < 1e-6, "flank_ratio_median = median(5,10)=7.5")
    check(int(row[5]) == 2, "flank n_samples_used = 2")
    check(row[6] == 100, "flank_len_used = candidate length")

    am = M.antisense_metrics(meta, c2, c1, cols, bams, regmap, min_reads=1)
    arow = [r for r in am if r[0] == "a1"][0]
    # ratios S1=5/50=0.1 ; S2=10/100=0.1 -> median 0.1
    check(abs(float(arow[6]) - 0.1) < 1e-6, "antisense ratio_median")
    check(arow[8] == 2, "antisense n = 2")


def test_intronic_relation():
    print("intronic relation (class i sense/antisense split)")
    gene_loc = {"G1": ("1", 1000, 1500, "+"), "G2": ("1", 1000, 1500, "-")}
    check(C.intronic_relation("i", "+", "G1", gene_loc) == "i_sense",
          "same-strand host -> i_sense")
    check(C.intronic_relation("i", "+", "G2", gene_loc) == "i_antisense",
          "opposite-strand host -> i_antisense")
    check(C.intronic_relation("i", "+", "GX", gene_loc) == "i_sense",
          "unknown host -> conservative i_sense")
    check(C.intronic_relation("i", "+", None, gene_loc) == "i_sense",
          "no host -> conservative i_sense")
    check(C.intronic_relation("x", "+", "G2", gene_loc) is None,
          "non-i class -> None")


def test_regions_host_fallback():
    print("regions host resolution (cmp_ref -> tx2gene) + host strand")
    d = tempfile.mkdtemp()
    ref = os.path.join(d, "ref.gtf")
    frozen = os.path.join(d, "frozen.gtf")
    rec = os.path.join(d, "rec.tsv")
    with open(ref, "w") as fh:
        fh.write('1\tX\tgene\t1000\t1500\t.\t+\t.\tgene_id "G1"; gene_biotype "protein_coding";\n')
        fh.write('1\tX\ttranscript\t1000\t1500\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n')
        fh.write('1\tX\texon\t1000\t1100\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n')
        fh.write('1\tX\texon\t1400\t1500\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n')
    with open(frozen, "w") as fh:
        # class x with only cmp_ref (reference transcript) -> host via tx2gene
        fh.write('1\tS\ttranscript\t1000\t1100\t.\t-\t.\tgene_id "M2"; transcript_id "M2.1"; class_code "x"; cmp_ref "T1"; cmp_ref_gene "G1sym";\n')
        fh.write('1\tS\texon\t1000\t1100\t.\t-\t.\tgene_id "M2"; transcript_id "M2.1";\n')
        # class i inside a host intron, cmp_ref only
        fh.write('1\tS\ttranscript\t1200\t1300\t.\t+\t.\tgene_id "M3"; transcript_id "M3.1"; class_code "i"; cmp_ref "T1"; cmp_ref_gene "G1sym";\n')
        fh.write('1\tS\texon\t1200\t1300\t.\t+\t.\tgene_id "M3"; transcript_id "M3.1";\n')
    with open(rec, "w") as fh:
        fh.write("ID\tn_reps_expressed\tpasses_recurrence\tclass_code\n")
        fh.write("T1\t5\tTRUE\t=\n")

    class A:
        pass
    a = A()
    a.tier1_classes = ["u", "x", "i", "j", "o"]
    a.tier1a_classes = ["u", "x", "i"]
    a.tier1_parent_biotypes = ["lncRNA"]
    a.annotated_lncrna_biotypes = ["lncRNA"]
    a.flank_bp = 500
    a.min_flank_len = 50
    a.cap = 100
    a.seed = 1
    a.expr_min_reps = 3
    out = os.path.join(d, "out")
    reference = C.load_reference(ref, ["lncRNA"])
    frozen_g = C.load_gtf(frozen, {"transcript", "exon"})
    R.build(reference, frozen_g, rec, out, a)
    _, regions = C.read_tsv(os.path.join(out, "regions.tsv"))
    ids = {r["element_id"]: r for r in regions}
    check("cand_antisense:M2.1" in ids, "class x candidate found via cmp_ref_gene")
    check(ids.get("cand_antisense:M2.1", {}).get("strand") == "+",
          "antisense region placed on the host strand")
    # class i host-intron flank activates
    flank_saf = os.path.join(out, "flank.saf")
    with open(flank_saf) as fh:
        lines = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
    left = [l for l in lines if l[0] == "cand:M3.1|left"]
    check(bool(left) and left[0][2] == "1101" and left[0][3] == "1399",
          "class i host-intron flank used (not genomic flanks)")


def test_thresholds():
    print("thresholds")
    check(T.decide(2.5, None, 2.0, 1000, 200, 0.95, 0.90)[:2] == (2.5, "derived"), "derive")
    check(T.decide(2.5, 3.0, 2.0, 1000, 200, 0.95, 0.90)[:2] == (3.0, "override"), "override wins")
    v, o, s = T.decide(2.5, None, 2.0, 10, 200, 0.95, 0.90)
    check((v, o, s) == (2.0, "default", "thresholds_unreliable"), "too few controls -> default")
    v, o, s = T.decide(2.5, None, 2.0, 1000, 200, 0.50, 0.90)
    check(s == "thresholds_unreliable", "bad separation -> unreliable")
    rows = [{"element_id": "e1", "set": "pos_flank", "norm_group": "g", "flank_ratio_median": "4"},
            {"element_id": "e1", "set": "pos_flank", "norm_group": "h", "flank_ratio_median": "6"},
            {"element_id": "e2", "set": "pos_flank", "norm_group": "g", "flank_ratio_median": "2"}]
    em = T.element_medians(rows, "flank_ratio_median")
    check(em["e1"] == 5.0 and em["e2"] == 2.0, "element medians over groups")


def test_tier_final():
    print("tier_final")
    check(F.flank_pass([3, 4], 2.5) is True, "flank pass true")
    check(F.flank_pass([1, 4], 2.5) is False, "flank pass false")
    check(F.flank_pass([], 2.5) is None, "flank pass unevaluable")
    check(F.antisense_risk([0.1, 0.3], 0.15) is True, "antisense risk")
    check(F.antisense_risk([0.1], 0.05) is False, "antisense no risk")
    check(F.antisense_status([0.1]) == "evaluated", "antisense_status evaluated")
    check(F.antisense_status([]) == "not_evaluable", "antisense_status not_evaluable")
    t1a, t1, pb = {"u", "x", "i"}, {"u", "x", "i", "j", "o"}, {"lncRNA"}
    # monoexonic u, allow=False -> demoted and marked
    nt, se, fp, ar, st, fin = F.evaluate_tx("M1", "u", None, True,
                                            {"cand:M1": [3.0]}, {}, 2.5, 0.15,
                                            t1, t1a, pb, False)
    check(se is True and st is False and fin is False,
          "monoexonic u demoted by categorical rule")
    # monoexonic u, allow=True -> tier_static; flank does not gate
    _, se2, _, _, st2, fin2 = F.evaluate_tx("M1", "u", None, True,
                                            {"cand:M1": [1.0]}, {}, 2.5, 0.15,
                                            t1, t1a, pb, True)
    check(se2 is False and st2 is True and fin2 is True,
          "monoexonic u kept when allow_single_exon=True (flank ignored)")
    # multi-exonic x, no flank needed -> final
    _, _, _, _, st3, fin3 = F.evaluate_tx("M3", "x", None, False,
                                          {}, {"cand_antisense:M3": [0.5]}, 2.5, 0.15,
                                          t1, t1a, pb, False)
    check(st3 is True and fin3 is True, "multiexonic x -> final")
    # multi-exonic x with antisense risk: D1 -> NOT dropped (F1.4 informational)
    _, _, _, ar4, st4, fin4 = F.evaluate_tx("M4", "x", None, False,
                                            {}, {"cand_antisense:M4": [0.1]}, 2.5, 0.15,
                                            t1, t1a, pb, False)
    check(st4 is True and fin4 is True and ar4 is True,
          "multiexonic x antisense risk -> tier_final = tier_static (D1)")
    # j with coding parent -> not static
    nt5, _, _, _, st5, fin5 = F.evaluate_tx("M5", "j", "protein_coding", False,
                                            {}, {}, 2.5, 0.15, t1, t1a, pb, False)
    check(nt5 is None and st5 is False and fin5 is False, "j coding parent -> not tier1")
    # no flank threshold (F1.3 disabled): flank_pass stays NA, tier_final unaffected
    _, _, fp6, _, _, fin6 = F.evaluate_tx("M6", "x", None, False,
                                          {}, {}, None, 0.15, t1, t1a, pb, False)
    check(fp6 is None and fin6 is True, "flank disabled -> multi-exonic x still final")


def test_class_composition():
    print("class composition (F1.6, D2 units)")
    d = tempfile.mkdtemp()
    ref = os.path.join(d, "ref.gtf")
    frozen = os.path.join(d, "frozen.gtf")
    with open(ref, "w") as fh:
        fh.write('1\tX\tgene\t100\t200\t.\t+\t.\tgene_id "L1"; gene_biotype "lncRNA";\n')
        fh.write('1\tX\tgene\t300\t400\t.\t-\t.\tgene_id "P1"; gene_biotype "protein_coding";\n')
    with open(frozen, "w") as fh:
        # gene M1: two transcripts, class i (multi) + class x (multi) -> Tier1 set {i,x}
        fh.write('1\tS\ttranscript\t100\t200\t.\t+\t.\tgene_id "M1"; transcript_id "M1.1"; class_code "i";\n')
        fh.write('1\tS\texon\t100\t120\t.\t+\t.\tgene_id "M1"; transcript_id "M1.1";\n')
        fh.write('1\tS\texon\t180\t200\t.\t+\t.\tgene_id "M1"; transcript_id "M1.1";\n')
        fh.write('1\tS\ttranscript\t100\t200\t.\t+\t.\tgene_id "M1"; transcript_id "M1.2"; class_code "x";\n')
        fh.write('1\tS\texon\t100\t120\t.\t+\t.\tgene_id "M1"; transcript_id "M1.2";\n')
        fh.write('1\tS\texon\t180\t200\t.\t+\t.\tgene_id "M1"; transcript_id "M1.2";\n')
        # gene M2: single mono-exonic u -> demoted only when allow_single_exon=false
        fh.write('1\tS\ttranscript\t100\t110\t.\t+\t.\tgene_id "M2"; transcript_id "M2.1"; class_code "u";\n')
        fh.write('1\tS\texon\t100\t110\t.\t+\t.\tgene_id "M2"; transcript_id "M2.1";\n')
        # gene M3: j with protein_coding parent -> Tier 2 always
        fh.write('1\tS\ttranscript\t300\t400\t.\t-\t.\tgene_id "M3"; transcript_id "M3.1"; class_code "j"; cmp_ref "P1"; cmp_ref_gene "P1sym";\n')
        fh.write('1\tS\texon\t300\t320\t.\t-\t.\tgene_id "M3"; transcript_id "M3.1";\n')
        fh.write('1\tS\texon\t380\t400\t.\t-\t.\tgene_id "M3"; transcript_id "M3.1";\n')
    t1a, t1, pb = {"u", "x", "i"}, {"u", "x", "i", "j", "o"}, {"lncRNA"}
    gene_biotype = {"L1": "lncRNA", "P1": "protein_coding"}
    frozen_g = C.load_gtf(frozen, {"transcript", "exon"})
    exons = C.exon_map(frozen_g)
    tx_off, genes_off = F.class_composition(frozen_g, exons, gene_biotype, t1, t1a, pb, False)
    check(tx_off == {"i": 1, "x": 1}, "transcripts Tier1 by class (allow=false)")
    check(sum(genes_off.values()) == 1, "M2 demoted: only M1 keeps Tier1 (allow=false)")
    check(genes_off[frozenset({"i", "x"})] == 1, "gene class set {i,x} counted once")
    tx_on, genes_on = F.class_composition(frozen_g, exons, gene_biotype, t1, t1a, pb, True)
    check(tx_on == {"i": 1, "x": 1, "u": 1}, "mono u added when allow=true")
    check(sum(genes_on.values()) == 2, "M2 recovered when allow=true")
    # j with coding parent never enters the composition
    check(not any("j" in cs for cs in genes_on), "j coding parent excluded")


def test_merge_and_assignment():
    print("merge and assignment")
    d = tempfile.mkdtemp()
    def one(path, bam, val):
        with open(path, "w") as fh:
            fh.write("# Program:featureCounts\n")
            fh.write("Geneid\tChr\tStart\tEnd\tStrand\tLength\t%s\n" % bam)
            fh.write("e1|candidate\t1\t100\t199\t+\t100\t%d\n" % val)
    one(os.path.join(d, "a.s2.tsv"), "bamA.bam", 10)
    one(os.path.join(d, "b.s2.tsv"), "bamB.bam", 20)
    meta, counts, cols = M.merge_featurecounts([os.path.join(d, "b.s2.tsv"),
                                                os.path.join(d, "a.s2.tsv")])
    check(cols == ["bamA.bam", "bamB.bam"], "merge sorts columns by bam")
    check(counts["e1|candidate"] == [10, 20], "merge aligns counts to sorted columns")
    bams = [("S1", "g1", "bamA.bam"), ("S2", "g1", "bamB.bam")]
    sg = M.check_assignment(cols, bams)
    check(sg == {"S1": "g1", "S2": "g1"}, "assignment resolves one group per sample")
    try:
        M.check_assignment(["unknown.bam"], bams)
        check(False, "unknown column must fail")
    except SystemExit:
        check(True, "unknown column fails")
    try:
        M.check_assignment(cols, [("S1", "g1", "bamA.bam"), ("S1", "g2", "bamA.bam")])
        check(False, "conflicting groups must fail")
    except SystemExit:
        check(True, "sample with two norm_groups fails")


def test_featurecounts_semicolon_fixture():
    print("featureCounts semicolon fixture")
    fx = os.path.join(HERE, "fixtures", "featurecounts", "anti_s2.tsv")
    meta, counts, cols = M.parse_featurecounts(fx)
    check(cols == ["NB4_Exp4_Scramble.markdup.sorted.bam"], "fixture has one sample column")
    gid = "pos_antisense:ENST00000359906|region"
    check(meta[gid]["Length"] == 736, "semicolon row Length parsed")
    check(meta[gid]["Start"] == "42509028;42510631;42511816;42512118",
          "semicolon Start kept as text")
    check(meta[gid]["Strand"] == "+;+;+;+", "semicolon Strand kept as text")
    check(counts[gid] == [71], "semicolon row count parsed")
    gid0 = "pos_antisense:ENST00000316807|region"
    check(meta[gid0]["Length"] == 892 and counts[gid0] == [0], "zero-count row parsed")
    _, s1_counts, _ = M.parse_featurecounts(
        os.path.join(HERE, "fixtures", "featurecounts", "anti_s1.tsv"))
    check(s1_counts[gid] == [634], "s1 fixture parses the same gid")


def _fc1(path, bam, rows):
    """Single-sample featureCounts table (one count column)."""
    with open(path, "w") as fh:
        fh.write("# Program:featureCounts\n")
        fh.write("Geneid\tChr\tStart\tEnd\tStrand\tLength\t%s\n" % bam)
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


def test_merge_multi_sample_duplicate():
    print("merge multi-sample + duplicate GeneID")
    d = tempfile.mkdtemp()
    a = os.path.join(d, "a.s2.tsv")
    b = os.path.join(d, "b.s2.tsv")
    _fc1(a, "bamA.bam", [
        ("g1|region", "1;1", "10;20", "15;25", "+;+", 30, 4),
        ("g1|region", "1", "10", "25", "+", 20, 6),
        ("g2|region", "1", "50", "80", "-", 31, 9),
    ])
    _fc1(b, "bamB.bam", [
        ("g1|region", "1", "10", "25", "+", 50, 7),
        ("g2|region", "1", "50", "80", "-", 31, 3),
    ])
    meta, counts, cols = M.merge_featurecounts([b, a])
    check(cols == ["bamA.bam", "bamB.bam"], "merge sorts by bam name")
    check(counts["g1|region"] == [10, 7], "duplicate GeneID counts summed per sample")
    check(meta["g1|region"]["Length"] == 50, "duplicate GeneID lengths summed")
    check(counts["g2|region"] == [9, 3], "second feature carried through")


def test_cardinality_guard():
    print("cardinality guard")
    d = tempfile.mkdtemp()
    s2 = os.path.join(d, "only.s2.tsv")
    s1 = os.path.join(d, "only.s1.tsv")
    _fc1(s2, "only.bam", [("a1|region", "1", "100", "199", "+", 100, 50)])
    _fc1(s1, "only.bam", [("a1|region", "1", "100", "199", "+", 100, 5)])
    regions = os.path.join(d, "regions.tsv")
    C.write_tsv([("a1", "antisense", "candidate", "1", "+", "")],
                regions, ["element_id", "kind", "set", "chrom", "strand", "note"])
    coldata = os.path.join(d, "coldata.tsv")
    with open(coldata, "w") as fh:
        fh.write("S1\tgrp1\tonly.bam\nS2\tgrp1\tother.bam\n")
    out = os.path.join(d, "out.tsv")
    argv = ["f1_metrics.py", "antisense", "--s2", s2, "--s1", s1,
            "--regions", regions, "--coldata", coldata, "--out", out]
    old = sys.argv
    sys.argv = argv
    try:
        M.main()
        check(False, "partial cohort must fail")
    except SystemExit as e:
        check("incomplete featureCounts input" in str(e), "partial cohort fails clearly")
    finally:
        sys.argv = old


def main():
    for fn in (test_intervals, test_stats, test_attrs, test_regions_pure,
               test_intronic_relation,
               test_regions_integration, test_regions_host_fallback, test_metrics,
               test_merge_and_assignment,
               test_featurecounts_semicolon_fixture, test_merge_multi_sample_duplicate,
               test_cardinality_guard, test_thresholds, test_tier_final,
               test_class_composition):
        fn()
    print()
    if FAILS:
        print("FAILED: %d" % len(FAILS))
        for f in FAILS:
            print("  - %s" % f)
        sys.exit(1)
    print("ALL F1 LOGIC TESTS PASSED")


if __name__ == "__main__":
    main()
