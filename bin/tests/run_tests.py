#!/usr/bin/env python3
"""Run the bin/tests suite without Nextflow (F0.12a).

Executes every `test_*.py` in this directory in a separate Python process and, when
given the two result tables, also runs the post-run consistency check. Standard
library only; no R, no samtools, no Nextflow.

  python3 bin/tests/run_tests.py
  python3 bin/tests/run_tests.py \
      results/discovery/gene_catalog/gene_catalog.tsv \
      results/f1/tier_final/tier_final_genes.tsv

Fixtures live in `bin/tests/fixtures/`. The Nextflow end-to-end dataset
(`tests/data/`, `conf/test.config`) is provided by F0.12b.
"""

import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))


def run(cmd):
    print("== %s" % " ".join(cmd))
    proc = subprocess.run(cmd)
    return proc.returncode


def check_bin_executables():
    """Nextflow puts bin/ on PATH and runs the file directly, so every script
    invoked as a bare command must carry the executable bit. A missing +x fails
    at runtime with exit 126, which is invisible to the unit tests."""
    bin_dir = os.path.join(ROOT, "bin")
    nf_files = [os.path.join(ROOT, "main.nf")]
    for d in ("modules/local", "workflows"):
        full = os.path.join(ROOT, d)
        if os.path.isdir(full):
            nf_files += [os.path.join(full, f) for f in os.listdir(full)
                         if f.endswith(".nf")]

    invoked = set()
    for path in nf_files:
        if not os.path.isfile(path):
            continue
        with open(path) as fh:
            for line in fh:
                m = re.match(r"\s*([A-Za-z0-9_]+\.py)(\s|$)", line)
                if m:
                    invoked.add(m.group(1))

    bad = [n for n in sorted(invoked)
           if os.path.isfile(os.path.join(bin_dir, n))
           and not os.access(os.path.join(bin_dir, n), os.X_OK)]
    print("== bin/*.py executable check (%d invoked as commands)" % len(invoked))
    if bad:
        print("  FAIL not executable: %s" % ", ".join(bad))
        return 1
    print("  ok  all invoked scripts are executable")
    return 0


def main():
    rc = check_bin_executables()
    tests = sorted(f for f in os.listdir(HERE)
                   if f.startswith("test_") and f.endswith(".py"))
    if not tests:
        print("no test_*.py found in %s" % HERE)
        return 1
    for t in tests:
        rc |= run([sys.executable, os.path.join(HERE, t)])

    if len(sys.argv) == 3:
        rc |= run([sys.executable, os.path.join(HERE, "check_tier_consistency.py"),
                   sys.argv[1], sys.argv[2]])
    else:
        print("(skipping check_tier_consistency.py: pass gene_catalog.tsv and "
              "tier_final_genes.tsv to enable it)")

    print()
    print("SUITE %s" % ("PASSED" if rc == 0 else "FAILED"))
    return rc


if __name__ == "__main__":
    sys.exit(main())
