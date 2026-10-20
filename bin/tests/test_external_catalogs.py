#!/usr/bin/env python3
"""F2.2 external-catalog registration tests. Stdlib only.

Run: python3 bin/tests/test_external_catalogs.py
"""

import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import external_catalogs as E  # noqa: E402

FAILS = []


def check(cond, msg):
    print(("  ok  " if cond else "  FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def test_register():
    print("register")
    d = tempfile.mkdtemp()
    with open(os.path.join(d, "lncipedia.bed"), "w") as fh:
        fh.write("1\t99\t200\ta\t0\t+\n1\t300\t400\tb\t0\t-\n1\t500\t600\tc\t0\t+\n")
    with open(os.path.join(d, "noncode.bed"), "w") as fh:
        fh.write("1\t99\t200\td\t0\t+\n1\t300\t400\te\t0\t-\n")
    entries = [
        {"name": "lncipedia", "version": "5.2", "url": "http://example/lncipedia.bed", "file": os.path.join(d, "lncipedia.bed")},
        {"name": "noncode", "version": "v6", "url": "", "file": os.path.join(d, "noncode.bed")},
        {"name": "fantomcat", "version": "", "url": "", "file": os.path.join(d, "missing.bed")},
    ]
    rows = {r["name"]: r for r in E.register(entries, d)}
    check(rows["lncipedia"]["n_features"] == 3 and rows["lncipedia"]["present"] == "TRUE",
          "lncipedia 3 features present")
    check(len(rows["lncipedia"]["md5"]) == 32, "md5 computed")
    check(rows["noncode"]["n_features"] == 2, "noncode 2 features")
    check(rows["fantomcat"]["present"] == "FALSE" and rows["fantomcat"]["n_features"] == 0,
          "missing file -> present FALSE, 0 features")


def main():
    test_register()
    print()
    if FAILS:
        print("FAILED: %d" % len(FAILS))
        sys.exit(1)
    print("ALL EXTERNAL CATALOG TESTS PASSED")


if __name__ == "__main__":
    main()
