#!/usr/bin/env python3
"""F0.7 guard tests on synthetic RSeQC outputs. Stdlib only.

Run: python3 bin/tests/test_strandedness.py
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import check_strandedness as S  # noqa: E402

FIX = os.path.join(HERE, "fixtures", "rseqc")
FAILS = []


def check(cond, msg):
    print(("  ok  " if cond else "  FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def rec(name):
    with open(os.path.join(FIX, name)) as fh:
        return (name, S.parse_one(fh.read()))


def test_parse():
    print("parse")
    p = S.parse_one(open(os.path.join(FIX, "reverse_a.strandedness.txt")).read())
    check(p["data_type"] == "paired" and abs(p["reverse"] - 0.88) < 1e-9, "parse PE reverse")


def test_ok():
    print("reverse matches declared")
    errs, warns, rows = S.evaluate([rec("reverse_a.strandedness.txt"),
                                    rec("reverse_b.strandedness.txt")], "reverse", 0.70, False)
    check(errs == [], "reverse data with declared reverse passes")
    check(len(rows) == 2, "summary rows emitted")


def test_declared_forward():
    print("declared forward on reverse data")
    errs, _, _ = S.evaluate([rec("reverse_a.strandedness.txt")], "forward", 0.70, False)
    check(any("majority" in e for e in errs), "declared forward on reverse fails")


def test_declared_unstranded():
    print("declared unstranded on reverse data")
    errs, _, _ = S.evaluate([rec("reverse_a.strandedness.txt")], "unstranded", 0.70, False)
    check(any("unstranded" in e for e in errs), "declared unstranded on stranded data fails")


def test_pe_se_incoherent():
    print("mixed PE/SE across samples")
    errs, _, _ = S.evaluate([rec("reverse_a.strandedness.txt"),
                             rec("single_reverse.strandedness.txt")], "reverse", 0.70, False)
    check(any("mixed library types" in e for e in errs), "PE+SE mixture fails")


def test_low_fidelity():
    print("low fidelity warn vs strict")
    errs, warns, _ = S.evaluate([rec("low_fidelity.strandedness.txt")], "reverse", 0.70, False)
    check(errs == [] and any("fidelity floor" in w for w in warns), "low fidelity warns, no fail")
    errs2, _, _ = S.evaluate([rec("low_fidelity.strandedness.txt")], "reverse", 0.70, True)
    check(any("fidelity floor" in e for e in errs2), "low fidelity fails under strict")


def main():
    for fn in (test_parse, test_ok, test_declared_forward, test_declared_unstranded,
               test_pe_se_incoherent, test_low_fidelity):
        fn()
    print()
    if FAILS:
        print("FAILED: %d" % len(FAILS))
        sys.exit(1)
    print("ALL STRANDEDNESS GUARD TESTS PASSED")


if __name__ == "__main__":
    main()
