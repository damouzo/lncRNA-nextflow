#!/usr/bin/env python3
"""H11 / 0.1: catalog -> universe funnel tests. Stdlib only.

The headline case is a catalog gene absent from the DE universe (e.g. the 27 bp
ENSG00000310590, below Salmon's k=31) must surface as step=not_quantified.
"""

import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
SCRIPT = os.path.join(ROOT, "bin", "build_funnel.py")

fails = 0


def check(cond, msg):
    global fails
    print(("  ok  " if cond else "  FAIL ") + msg)
    if not cond:
        fails += 1


def write(path, text):
    with open(path, "w") as fh:
        fh.write(text)


d = tempfile.mkdtemp("funnel_")
# G_short stands for a <31 nt transcript that Salmon's k=31 index cannot quantify.
write(os.path.join(d, "gene_catalog.tsv"),
      "gene_id\torigin\ttier_static\n"
      "G_ok\tnovel\tTRUE\n"
      "G_short\tnovel\tTRUE\n"
      "G_ann\tannotated_lncrna\tNA\n")
# G_short is absent from the universe; G_ann is present but filtered out.
write(os.path.join(d, "MNC_universe.tsv"),
      "gene_id\ttotal_count\tkeep\tfilter_reason\n"
      "G_ok\t500\tTRUE\t\n"
      "G_ann\t2\tFALSE\tlow_expression\n")
write(os.path.join(d, "coding_consensus.tsv"),
      "ID\tmRNA\tCPAT_noncoding\tCPC2_noncoding\tconsensus_noncoding\n"
      "T1\tno\tTRUE\tTRUE\tTRUE\nT2\tyes\tFALSE\tFALSE\tFALSE\n")
write(os.path.join(d, "recurrence.tsv"),
      "transcript_id\tn_reps_expressed\tpasses_recurrence\n"
      "T1\t3\tTRUE\nT2\t1\tFALSE\n")

out = os.path.join(d, "out")
proc = subprocess.run([sys.executable, SCRIPT,
                       "--gene-catalog", os.path.join(d, "gene_catalog.tsv"),
                       "--universe", os.path.join(d, "MNC_universe.tsv"),
                       "--coding-consensus", os.path.join(d, "coding_consensus.tsv"),
                       "--recurrence", os.path.join(d, "recurrence.tsv"),
                       "--outdir", out], capture_output=True, text=True)
check(proc.returncode == 0, "build_funnel.py exits 0")

rows = []
with open(os.path.join(out, "funnel.tsv")) as fh:
    header = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))

check(header == ["step", "gene_id", "n", "reason"], "funnel.tsv has step/gene_id/n/reason")

agg = {r["step"]: r["n"] for r in rows if r["gene_id"] == ""}
check(agg.get("catalog") == "3", "catalog count")
check(agg.get("de_universe") == "2", "de_universe count")
check(agg.get("quantified_catalog") == "2", "quantified_catalog = catalog intersect universe")
check(agg.get("not_quantified") == "1", "one catalog gene absent from the universe")
check(agg.get("cpat_noncoding") == "1" and agg.get("coding_consensus_noncoding") == "1",
      "CPAT noncoding + consensus steps present")
check(agg.get("recurrence_pass") == "1", "recurrence step present")

short = [r for r in rows if r["gene_id"] == "G_short"]
check(len(short) == 1 and short[0]["step"] == "not_quantified"
      and short[0]["reason"] == "not_quantified",
      "the <31 nt gene comes out as not_quantified with a reason")

low = [r for r in rows if r["step"] == "universe_filtered_low_expression"]
check(len(low) == 1 and low[0]["n"] == "1", "low-expression filtered genes are counted")

print()
if fails:
    print("FAILED: %d" % fails)
    sys.exit(1)
print("ALL FUNNEL TESTS PASSED")
