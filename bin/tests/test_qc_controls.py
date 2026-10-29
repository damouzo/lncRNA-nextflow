#!/usr/bin/env python3
"""Synthetic tests for bin/qc_controls.py (positive-control report).

Cases: a control that goes down (pass), a control that goes up or is not
significant (fail), a primary-patient contrast (not_applicable), a control gene
missing from the DE table (fail). The report is diagnostic only: it never fails
the run, so every case exits 0.
"""

import csv
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
SCRIPT = os.path.join(ROOT, "bin", "qc_controls.py")

fails = 0


def check(cond, msg):
    global fails
    print(("  ok  " if cond else "  FAIL ") + msg)
    if not cond:
        fails += 1


def write_de(tmp, group, ct, gene_rows):
    """gene_rows: {gene_id: (lfc, padj)}."""
    name = "%s_DE.tsv" % ct
    with open(os.path.join(tmp, name), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene_id", "baseMean", "log2FoldChange", "lfcSE",
                    "pvalue", "padj", "contrast"])
        for gid, (lfc, padj) in gene_rows.items():
            w.writerow([gid, 100.0, lfc, 0.1, 0.001, padj, ct])
    return name


def run(args):
    return subprocess.run([sys.executable, SCRIPT] + args,
                          capture_output=True, text=True)


def read_out(tmp):
    rows = []
    with open(os.path.join(tmp, "qc_controls.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append(r)
    return {(r["group"], r["contrast"], r["control_gene"]): r for r in rows}


DDX41 = "ENSG00000183258"
DHX34 = "ENSG00000134815"

tmp = tempfile.mkdtemp("qcc_")

# contrasts: two shRNA contrasts + one primary-patient contrast
with open(os.path.join(tmp, "contrasts.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["NB4_DDX41sh1_vs_Scramble", "NB4_DDX41_sh1", "NB4_Scramble", ""])
    w.writerow(["NB4_DHX34sh1_vs_Scramble", "NB4_DHX34_sh1", "NB4_Scramble", ""])
    w.writerow(["MNC_DDX41_vs_Healthy", "MNC_DDX41_Patient", "MNC_Healthy", ""])

# pass: DDX41 down in its sh contrast
f1 = write_de(tmp, "NB4", "NB4_DDX41sh1_vs_Scramble", {DDX41: (-1.8, 0.001)})
# fail: DHX34 up in its sh contrast
f2 = write_de(tmp, "NB4", "NB4_DHX34sh1_vs_Scramble", {DHX34: (0.9, 0.001)})
# not_applicable: primary patient contrast, no sh
f3 = write_de(tmp, "MNC", "MNC_DDX41_vs_Healthy", {DDX41: (0.0, 0.9)})

with open(os.path.join(tmp, "manifest.tsv"), "w", newline="") as fh:
    fh.write("group\tcontrast\tfile\n")
    fh.write("NB4\tNB4_DDX41sh1_vs_Scramble\t%s\n" % f1)
    fh.write("NB4\tNB4_DHX34sh1_vs_Scramble\t%s\n" % f2)
    fh.write("MNC\tMNC_DDX41_vs_Healthy\t%s\n" % f3)

# a control failure must not fail the process, but must be recorded
res = run(["--manifest", os.path.join(tmp, "manifest.tsv"),
           "--de-dir", tmp, "--contrasts", os.path.join(tmp, "contrasts.tsv"),
           "--outdir", tmp])
check(res.returncode == 0, "a control failure exits 0 (report-only)")
out = read_out(tmp)
check(out[("NB4", "NB4_DDX41sh1_vs_Scramble", "DDX41")]["status"] == "pass",
      "DDX41 down -> pass")
check(out[("NB4", "NB4_DHX34sh1_vs_Scramble", "DHX34")]["status"] == "fail",
      "DHX34 up -> fail")
check(out[("NB4", "NB4_DDX41sh1_vs_Scramble", "DHX34")]["status"] == "not_applicable",
      "DHX34 in a DDX41 contrast -> not_applicable")
check(out[("MNC", "MNC_DDX41_vs_Healthy", "DDX41")]["status"] == "not_applicable",
      "primary-patient contrast -> not_applicable")
check("NOTE" in res.stderr, "a control failure is reported on stderr")

# missing control gene -> fail
tmp2 = tempfile.mkdtemp("qcc2_")
sub2 = os.path.join(tmp2, "NB4_deseq2_results")
os.makedirs(sub2)
f4 = write_de(sub2, "NB4", "NB4_DDX41sh1_vs_Scramble", {"ENSG00000000003": (-3.0, 0.001)})
with open(os.path.join(tmp2, "contrasts.tsv"), "w", newline="") as fh:
    fh.write("NB4_DDX41sh1_vs_Scramble\tNB4_DDX41_sh1\tNB4_Scramble\t\n")
res3 = run(["--de", os.path.join(sub2, f4),
            "--contrasts", os.path.join(tmp2, "contrasts.tsv"),
            "--outdir", tmp2])
out3 = read_out(tmp2)
check(out3[("NB4", "NB4_DDX41sh1_vs_Scramble", "DDX41")]["status"] == "fail",
      "control gene absent from the DE table -> fail")

print()
if fails:
    print("FAILED: %d" % fails)
    sys.exit(1)
print("ALL QC_CONTROLS TESTS PASSED")
