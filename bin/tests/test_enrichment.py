#!/usr/bin/env python3
"""1.4 tests: functional enrichment keeps full tables and tags the foreground.

clusterProfiler is not available offline, so this locks the module contract (ORA
only, no p/q cutoff, single foreground label) and checks that the per-contrast split
preserves the tagging column.
"""

import csv
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
BIN = os.path.dirname(HERE)

fails = 0


def check(cond, msg):
    global fails
    print(("  ok  " if cond else "  FAIL ") + msg)
    if not cond:
        fails += 1


src = open(os.path.join(ROOT, "modules", "local", "functional_enrichment.nf")).read()
check("pvalueCutoff = 1" in src and "qvalueCutoff = 1" in src,
      "ORA keeps every tested term (pvalueCutoff/qvalueCutoff = 1)")
check("gseGO" not in src and "GSEA" not in src and "gsea" not in src,
      "GSEA removed from the module")
check("de_lncrna_neighbours" in src and "significant_cis_targets" not in src,
      "single ORA foreground (de_lncrna_neighbours)")

# The per-contrast split must pass the tagging columns through.
d = tempfile.mkdtemp()
ora = os.path.join(d, "g_enrichment_ora.tsv")
with open(ora, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["ID", "Description", "p.adjust", "contrast", "input_set", "foreground_type"])
    w.writerow(["GO:1", "x", "0.5", "g__c1", "de_lncrna_neighbours", "de_lncrna_neighbours"])
    w.writerow(["GO:2", "y", "0.9", "g__c2", "de_lncrna_neighbours", "de_lncrna_neighbours"])

out = os.path.join(d, "c1")
rc = subprocess.run([sys.executable, os.path.join(BIN, "split_enrichment.py"),
                     "--ora", ora, "--contrast", "g__c1",
                     "--outdir", out]).returncode
check(rc == 0, "split_enrichment exits 0")
with open(os.path.join(out, "enrichment_ora.tsv")) as fh:
    header = fh.readline().rstrip("\n").split("\t")
    row = fh.readline().rstrip("\n").split("\t")
check("input_set" in header, "split keeps the input_set column")
check(row[header.index("input_set")] == "de_lncrna_neighbours",
      "split keeps the input_set value")

print()
if fails:
    print("FAILED: %d" % fails)
    sys.exit(1)
print("ALL ENRICHMENT TESTS PASSED")
