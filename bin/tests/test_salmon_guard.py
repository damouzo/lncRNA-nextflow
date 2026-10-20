#!/usr/bin/env python3
"""F0.9 guard tests on synthetic Salmon JSON fixtures. Stdlib only.

Run: python3 bin/tests/test_salmon_guard.py
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import check_salmon_libtype as G  # noqa: E402

FIX = os.path.join(HERE, "fixtures", "salmon")
FAILS = []


def check(cond, msg):
    print(("  ok  " if cond else "  FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def test_ok_isr():
    print("ok_isr (correct ISR)")
    errs, qc = G.check(os.path.join(FIX, "ok_isr"), declared="ISR")
    check(errs == [], "ISR correct passes (no hard error)")
    cf = qc["consistent_fraction"]
    check(cf is not None and abs(cf - 40906101 / 42358422) < 1e-9,
          "consistent_fraction reported as QC")
    check(abs(qc["compatible_fragment_ratio"] - 1.0) < 1e-9,
          "compatible_fragment_ratio reported as QC")


def test_wrong_libtype_isf():
    print("wrong_libtype_isf (declared ISF, contract ISR)")
    errs, _ = G.check(os.path.join(FIX, "wrong_libtype_isf"), declared="ISR")
    check(any(e.startswith("H4") for e in errs), "cmd ISF vs contract ISR fails (H4)")


def test_archived_iu():
    print("archived_iu (declared SR, expected IU, PE data)")
    errs, _ = G.check(os.path.join(FIX, "archived_iu"), declared="ISR")
    check(any(e.startswith("H1") for e in errs), "expected IU != declared SR fails (H1)")
    check(any(e.startswith("H3") for e in errs), "SR on paired data fails (H3)")


def test_pe_se_mismatch():
    print("pe_se_mismatch (single read file, paired libtype)")
    errs, _ = G.check(os.path.join(FIX, "pe_se_mismatch"), declared="ISR")
    check(any(e.startswith("H3") for e in errs), "paired libtype with 1 read file fails (H3)")


def main():
    for fn in (test_ok_isr, test_wrong_libtype_isf, test_archived_iu, test_pe_se_mismatch):
        fn()
    print()
    if FAILS:
        print("FAILED: %d" % len(FAILS))
        sys.exit(1)
    print("ALL SALMON GUARD TESTS PASSED")


if __name__ == "__main__":
    main()
