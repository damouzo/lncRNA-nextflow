#!/usr/bin/env python3
"""F6.1 / F6.2 (recovery_analysis) tests. Stdlib only.

Run: python3 bin/tests/test_recovery_analysis.py
"""

import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.dirname(HERE)
sys.path.insert(0, BIN)

import f1_common as C  # noqa: E402
import recovery_analysis as R  # noqa: E402

FAILS = []


def check(cond, msg):
    print(("  ok  " if cond else "  FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def write(path, text):
    with open(path, "w") as fh:
        fh.write(text)


def make_catalog(d):
    p = os.path.join(d, "catalog_genes.tsv")
    write(p, "gene_id\torigin\ttier_final\n"
             "R1\tannotated_lncrna\tFALSE\n"
             "R2\tannotated_lncrna\tFALSE\n"
             "N1\tnovel\tTRUE\n"
             "N2\tnovel\tFALSE\n"
             "C1\treference_other\tFALSE\n")
    return p


def make_tpm(d, name, cols, values):
    p = os.path.join(d, name)
    lines = ["gene_id\t" + "\t".join(cols)]
    for gid, vals in values.items():
        lines.append(gid + "\t" + "\t".join(str(v) for v in vals))
    write(p, "\n".join(lines) + "\n")
    return p


def fixtures(d):
    cat = make_catalog(d)
    # expr pooled over 4 samples: R1 3x>=1 (expressed), R2 1x, N1 4x, N2 0x.
    tpm_a = make_tpm(d, "gA_tpm.tsv", ["s1", "s2"],
                     {"R1": [5, 5], "R2": [1, 0], "N1": [5, 5],
                      "N2": [0, 0], "C1": [9, 9]})
    tpm_b = make_tpm(d, "gB_tpm.tsv", ["s3", "s4"],
                     {"R1": [5, 0], "R2": [0, 0], "N1": [5, 5],
                      "N2": [0, 0], "C1": [9, 9]})
    return cat, [tpm_a, tpm_b]


def test_recovery_and_comparison():
    d = tempfile.mkdtemp()
    cat, tpms = fixtures(d)
    _, genes = C.read_tsv(cat)
    expr = R.load_expression(tpms)

    rec = dict(R.build_recovery(genes, expr, 1.0, 3))
    check(rec["reference_lncrna_genes"] == "2", "F6.1 counts reference lncRNAs")
    check(rec["reference_expressed"] == "1", "F6.1 pools samples across TPM tables")
    check(rec["reference_fraction_recovered"] == "0.5000", "F6.1 fraction recovered")
    check(rec["n_samples_pooled"] == "4", "F6.1 pools 4 samples")

    rows = {r[0]: r for r in R.build_comparison(genes, expr, 1.0, 3)}
    check(rows["novel_tier1"][1] == 1 and rows["novel_tier1"][2] == 1,
          "F6.2 novel_tier1 counted and expressed")
    check(rows["novel_tier2"][1] == 1 and rows["novel_tier2"][2] == 0,
          "F6.2 novel_tier2 counted, not expressed")
    check(rows["annotated_lncrna"][1] == 2 and rows["annotated_lncrna"][2] == 1,
          "F6.2 annotated lncRNAs counted and expressed")
    check(len(rows) == 3, "F6.2 ignores reference_other")


def test_min_samples_gate():
    d = tempfile.mkdtemp()
    cat, tpms = fixtures(d)
    _, genes = C.read_tsv(cat)
    expr = R.load_expression(tpms)
    rec = dict(R.build_recovery(genes, expr, 1.0, 4))
    check(rec["reference_expressed"] == "0",
          "min_samples gate excludes genes below the threshold")


def test_cli():
    d = tempfile.mkdtemp()
    cat, tpms = fixtures(d)
    out = os.path.join(d, "out")
    cmd = [sys.executable, os.path.join(BIN, "recovery_analysis.py"),
           "--catalog", cat, "--outdir", out]
    for t in tpms:
        cmd += ["--tpm", t]
    rc = subprocess.run(cmd).returncode
    check(rc == 0, "CLI exits 0")
    for f in ("recovery_reference.tsv", "tier1_vs_annotated.tsv"):
        check(os.path.exists(os.path.join(out, f)), "CLI wrote %s" % f)


def main():
    for t in (test_recovery_and_comparison, test_min_samples_gate, test_cli):
        t()
    print()
    if FAILS:
        print("RECOVERY TESTS FAILED: %d" % len(FAILS))
        return 1
    print("ALL RECOVERY TESTS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
