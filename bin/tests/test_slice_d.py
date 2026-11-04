#!/usr/bin/env python3
"""Slice D tests: split_contrasts (several contrasts per group) and the F2.3
novelty guard. Stdlib only, no Nextflow.

Run: python3 bin/tests/test_slice_d.py
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


def run(script, *args):
    return subprocess.run([sys.executable, os.path.join(BIN, script)] + list(args),
                          capture_output=True, text=True)


def test_split_contrasts():
    d = tempfile.mkdtemp()
    write(os.path.join(d, "coldata_MNC.csv"),
          "sample,condition,batch,norm_group\nS1,A,b1,MNC\nS2,B,b1,MNC\n")
    write(os.path.join(d, "coldata_MSC.csv"),
          "sample,condition,batch,norm_group\nS3,A,b2,MSC\nS4,B,b2,MSC\nS5,C,b2,MSC\n")
    write(os.path.join(d, "contrasts.tsv"),
          "A_vs_B\tA\tB\t\nA_vs_C\tA\tC\t\nA_vs_Z\tA\tZ\t\n")

    out = os.path.join(d, "out", "MNC")
    p = run("split_contrasts.py", "--coldata", os.path.join(d, "coldata_MNC.csv"),
            "--contrasts", os.path.join(d, "contrasts.tsv"), "--group", "MNC",
            "--outdir", out)
    check(p.returncode == 0, "split_contrasts MNC runs")
    files = sorted(os.listdir(out))
    check(files == ["MNC__A_vs_B.tsv"], "MNC keeps only the applicable contrast (%s)" % files)

    out2 = os.path.join(d, "out", "MSC")
    p = run("split_contrasts.py", "--coldata", os.path.join(d, "coldata_MSC.csv"),
            "--contrasts", os.path.join(d, "contrasts.tsv"), "--group", "MSC",
            "--outdir", out2)
    files2 = sorted(os.listdir(out2))
    check(p.returncode == 0 and files2 == ["MSC__A_vs_B.tsv", "MSC__A_vs_C.tsv"],
          "MSC keeps its two applicable contrasts (%s)" % files2)
    with open(os.path.join(out2, "MSC__A_vs_B.tsv")) as fh:
        body = fh.read()
    check(body.splitlines()[0] == "contrast_name\tnumerator\tdenominator\tbatch",
          "single-contrast file carries the header")

    write(os.path.join(d, "bad.tsv"), "A__B\tA\tB\t\n")
    p = run("split_contrasts.py", "--coldata", os.path.join(d, "coldata_MNC.csv"),
            "--contrasts", os.path.join(d, "bad.tsv"), "--group", "MNC",
            "--outdir", os.path.join(d, "out", "bad"))
    check(p.returncode != 0 and "__" in p.stderr, "contrast name with '__' is rejected")

    write(os.path.join(d, "slash.tsv"), "A/B\tA\tB\t\n")
    p = run("split_contrasts.py", "--coldata", os.path.join(d, "coldata_MNC.csv"),
            "--contrasts", os.path.join(d, "slash.tsv"), "--group", "MNC",
            "--outdir", os.path.join(d, "out", "slash"))
    check(p.returncode != 0 and "/" in p.stderr, "contrast name with '/' is rejected")

    # Names are interpolated into downstream shell commands: keep them quotable-free
    write(os.path.join(d, "space.tsv"), "A vs B\tA\tB\t\n")
    p = run("split_contrasts.py", "--coldata", os.path.join(d, "coldata_MNC.csv"),
            "--contrasts", os.path.join(d, "space.tsv"), "--group", "MNC",
            "--outdir", os.path.join(d, "out", "space"))
    check(p.returncode != 0 and "A vs B" in p.stderr,
          "contrast name with a space is rejected")

    write(os.path.join(d, "subst.tsv"), "A$(id)\tA\tB\t\n")
    p = run("split_contrasts.py", "--coldata", os.path.join(d, "coldata_MNC.csv"),
            "--contrasts", os.path.join(d, "subst.tsv"), "--group", "MNC",
            "--outdir", os.path.join(d, "out", "subst"))
    check(p.returncode != 0 and "A$(id)" in p.stderr,
          "contrast name with a command substitution is rejected")

    p = run("split_contrasts.py", "--coldata", os.path.join(d, "coldata_MSC.csv"),
            "--contrasts", os.path.join(d, "contrasts.tsv"), "--group", "MSC line",
            "--outdir", os.path.join(d, "out", "badgroup"))
    check(p.returncode != 0 and "group" in p.stderr,
          "unsafe group name is rejected")


def test_novelty():
    d = tempfile.mkdtemp()
    write(os.path.join(d, "consensus.gtf"),
          'chr1\ts\ttranscript\t1\t100\t.\t+\t.\ttranscript_id "T1";\n')
    write(os.path.join(d, "overlap.tsv"),
          "ID\tknown_catalog_hit\nT1\tTRUE\nMSTRG.1\tFALSE\n")

    ok = os.path.join(d, "ok.bed")
    write(ok, "chr1\t0\t100\tx\nchr1\t200\t300\ty\n")

    # no catalogs -> not_evaluated, no failure
    out = os.path.join(d, "n0.tsv")
    p = run("check_novelty.py", "--overlap", os.path.join(d, "overlap.tsv"),
            "--consensus-gtf", os.path.join(d, "consensus.gtf"), "--out", out)
    check(p.returncode == 0, "no catalogs does not fail")
    with open(out) as fh:
        check("not_evaluated" in fh.read(), "no catalogs -> novelty_status=not_evaluated")

    # catalogs present, hits -> evaluated
    out1 = os.path.join(d, "n1.tsv")
    p = run("check_novelty.py", "--overlap", os.path.join(d, "overlap.tsv"),
            "--consensus-gtf", os.path.join(d, "consensus.gtf"),
            "--bed", "lncipedia=" + ok, "--out", out1)
    check(p.returncode == 0, "catalogs with hits -> passes")
    with open(out1) as fh:
        check("novelty_status\tnovel_for_all_catalogs" in fh.read(),
              "annotated-only hits -> novel_for_all_catalogs")

    # chr-prefix alias: a BED without the prefix still matches the chr-prefixed consensus
    noprefix = os.path.join(d, "noprefix.bed")
    write(noprefix, "1\t0\t100\tx\nchrM\t0\t50\tm\n")
    out5 = os.path.join(d, "n5.tsv")
    p = run("check_novelty.py", "--overlap", os.path.join(d, "overlap.tsv"),
            "--consensus-gtf", os.path.join(d, "consensus.gtf"),
            "--bed", "noncode=" + noprefix, "--out", out5)
    check(p.returncode == 0, "non-prefixed BED matches chr-prefixed consensus (chr/M alias)")
    with open(out5) as fh:
        check("novelty_status\tnovel_for_all_catalogs" in fh.read(),
              "alias-normalized BED -> novel_for_all_catalogs")

    # a novel (MSTRG) candidate that hits a catalog -> known
    known = os.path.join(d, "known.tsv")
    write(known, "ID\tknown_catalog_hit\nT1\tTRUE\nMSTRG.1\tTRUE\n")
    out6 = os.path.join(d, "n6.tsv")
    p = run("check_novelty.py", "--overlap", known,
            "--consensus-gtf", os.path.join(d, "consensus.gtf"),
            "--bed", "lncipedia=" + ok, "--out", out6)
    check(p.returncode == 0, "novel candidate with a catalog hit -> passes")
    with open(out6) as fh:
        check("novelty_status\tknown" in fh.read(), "novel hit -> novelty_status=known")

    # empty bed -> fail
    empty = os.path.join(d, "empty.bed")
    write(empty, "")
    p = run("check_novelty.py", "--overlap", os.path.join(d, "overlap.tsv"),
            "--consensus-gtf", os.path.join(d, "consensus.gtf"),
            "--bed", "lncipedia=" + empty, "--out", os.path.join(d, "n2.tsv"))
    check(p.returncode != 0 and "0 features" in p.stderr, "BED with 0 features fails")

    # no chromosome in common -> fail
    other = os.path.join(d, "other.bed")
    write(other, "chrZ\t0\t100\tz\n")
    p = run("check_novelty.py", "--overlap", os.path.join(d, "overlap.tsv"),
            "--consensus-gtf", os.path.join(d, "consensus.gtf"),
            "--bed", "lncipedia=" + other, "--out", os.path.join(d, "n3.tsv"))
    check(p.returncode != 0 and "no chromosome" in p.stderr,
          "BED with no shared chromosome fails")

    # catalogs present but zero hit -> fail
    write(os.path.join(d, "nohit.tsv"), "ID\tknown_catalog_hit\nT1\tFALSE\n")
    p = run("check_novelty.py", "--overlap", os.path.join(d, "nohit.tsv"),
            "--consensus-gtf", os.path.join(d, "consensus.gtf"),
            "--bed", "lncipedia=" + ok, "--out", os.path.join(d, "n4.tsv"))
    check(p.returncode != 0 and "no known_catalog_hit" in p.stderr,
          "catalogs with zero hits fail")


def main():
    for fn in (test_split_contrasts, test_novelty):
        fn()
    print()
    if FAILS:
        print("FAILED: %d" % len(FAILS))
        for f in FAILS:
            print("  - %s" % f)
        sys.exit(1)
    print("ALL SLICE D TESTS PASSED")


if __name__ == "__main__":
    main()
