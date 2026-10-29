#!/usr/bin/env python3
"""1.4 tests: functional enrichment keeps full tables and tags the foreground.

clusterProfiler is not available offline, so this locks the module contract
(no p/q cutoff, ranking by the unshrunken `stat`, `foreground_type` column) and
checks that the per-contrast split preserves the new column.
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
check("pvalueCutoff = 1" in src.split("gseGO")[1] if "gseGO" in src else False,
      "GSEA keeps every tested term")
check('stat_col <- if ("stat" %in% colnames(pc_data))' in src,
      "GSEA ranks by the unshrunken `stat` (fallback to log2FoldChange)")
check("foreground_type" in src and "all_cis_targets_expressed" in src,
      "ORA tags the foreground (significant_cis_targets / all_cis_targets_expressed)")

# The per-contrast split must pass the new column through.
d = tempfile.mkdtemp()
ora = os.path.join(d, "g_enrichment_ora.tsv")
gsea = os.path.join(d, "g_enrichment_gsea.tsv")
with open(ora, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["ID", "Description", "p.adjust", "contrast", "foreground_type"])
    w.writerow(["GO:1", "x", "0.5", "g__c1", "significant_cis_targets"])
    w.writerow(["GO:2", "y", "0.9", "g__c2", "all_cis_targets_expressed"])
with open(gsea, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["ID", "NES", "contrast"])
    w.writerow(["GO:1", "1.5", "g__c1"])

out = os.path.join(d, "c1")
rc = subprocess.run([sys.executable, os.path.join(BIN, "split_enrichment.py"),
                     "--ora", ora, "--gsea", gsea, "--contrast", "g__c1",
                     "--outdir", out]).returncode
check(rc == 0, "split_enrichment exits 0")
with open(os.path.join(out, "enrichment_ora.tsv")) as fh:
    header = fh.readline().rstrip("\n").split("\t")
    row = fh.readline().rstrip("\n").split("\t")
check("foreground_type" in header, "split keeps the foreground_type column")
check(row[header.index("foreground_type")] == "significant_cis_targets",
      "split keeps the foreground_type value")

print()
if fails:
    print("FAILED: %d" % fails)
    sys.exit(1)
print("ALL ENRICHMENT TESTS PASSED")
