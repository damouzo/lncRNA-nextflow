#!/usr/bin/env python3
"""derive_seq_batch.py: FASTQ paths -> a canonical per-sample seq_batch value.

Covers the three real shapes in the cohort: a two-lane sample split across two
flowcells (the common case), a single-lane sample (A16 has no L1), and the
suffix-R replacement samples.
"""

import csv
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
SCRIPT = os.path.join(ROOT, "bin", "derive_seq_batch.py")

fails = 0


def check(cond, msg):
    global fails
    print(("  ok  " if cond else "  FAIL ") + msg)
    if not cond:
        fails += 1


d = tempfile.mkdtemp("seqbatch_")
sheet = os.path.join(d, "samplesheet.csv")
with open(sheet, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["sample", "condition", "batch", "norm_group",
                "fastq_1", "fastq_2", "bam", "seq_batch"])
    # Two lanes, two flowcells (the common case).
    w.writerow(["S_multi", "A", "Exp1", "G",
                "/x/A1/A1_LIB_23HJT3LT4_L8_1.fq.gz", "/x/A1/A1_LIB_23HJT3LT4_L8_2.fq.gz", "/x/A1.bam", "main"])
    w.writerow(["S_multi", "A", "Exp1", "G",
                "/x/A1/A1_LIB_23HJVTLT4_L1_1.fq.gz", "/x/A1/A1_LIB_23HJVTLT4_L1_2.fq.gz", "/x/A1.bam", "main"])
    # Single lane only.
    w.writerow(["S_single", "B", "Exp2", "G",
                "/x/A16/A16_LIB_23HJT3LT4_L8_1.fq.gz", "/x/A16/A16_LIB_23HJT3LT4_L8_2.fq.gz", "/x/A16.bam", "L4"])
    # Replacement sample: suffix R and a second flowcell.
    w.writerow(["S_rep", "C", "Exp3", "G",
                "/x/A39R/A39R_LIB_23HK3GLT4_L5_1.fq.gz", "/x/A39R/A39R_LIB_23HK3GLT4_L5_2.fq.gz", "/x/A39R.bam", "L5"])

subprocess.run([sys.executable, SCRIPT, "--samplesheet", sheet], check=True)

with open(sheet, newline="") as fh:
    rows = {r["sample"]: r["seq_batch"] for r in csv.DictReader(fh)}

check(rows.get("S_multi") == "23HJT3LT4_L8;23HJVTLT4_L1",
      "two-lane sample -> sorted, unique flowcell_lane set (%s)" % rows.get("S_multi"))
check(rows.get("S_single") == "23HJT3LT4_L8",
      "single-lane sample -> its only lane (%s)" % rows.get("S_single"))
check(rows.get("S_rep") == "23HK3GLT4_L5",
      "suffix-R replacement sample parsed (%s)" % rows.get("S_rep"))

# Idempotent: a second run must not change anything.
subprocess.run([sys.executable, SCRIPT, "--samplesheet", sheet], check=True)
with open(sheet, newline="") as fh:
    rows2 = {r["sample"]: r["seq_batch"] for r in csv.DictReader(fh)}
check(rows == rows2, "re-running is idempotent")

# The column is inserted when missing, right after batch.
sheet2 = os.path.join(d, "no_seq_batch.csv")
with open(sheet2, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["sample", "condition", "batch", "fastq_1", "fastq_2", "bam"])
    w.writerow(["S1", "A", "Exp1", "/x/A2/A2_LIB_23HJT3LT4_L8_1.fq.gz",
                "/x/A2/A2_LIB_23HJT3LT4_L8_2.fq.gz", "/x/A2.bam"])
subprocess.run([sys.executable, SCRIPT, "--samplesheet", sheet2], check=True)
with open(sheet2, newline="") as fh:
    hdr = next(csv.reader(fh))
check(hdr == ["sample", "condition", "batch", "seq_batch", "fastq_1", "fastq_2", "bam"],
      "seq_batch inserted after batch when absent (%s)" % hdr)

print()
if fails:
    print("FAILED: %d" % fails)
    sys.exit(1)
print("ALL DERIVE_SEQ_BATCH TESTS PASSED")
