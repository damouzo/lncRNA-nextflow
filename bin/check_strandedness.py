#!/usr/bin/env python3
"""F0.7 guard: empirical strandedness vs the `library_strandedness` contract.

Pure Python (stdlib), unit-tested with synthetic RSeQC `infer_experiment` outputs.
This is the orientation guard: `check_salmon_libtype.py` only echoes the `-l`
parameter, so a declaration that contradicts the data is caught here.

Fails when:
  - a sample's orientation cannot be parsed;
  - samples disagree on paired/single-end;
  - the empirical majority does not match `library_strandedness` (declaring
    forward / unstranded on reverse data stops the run);
  - `strict` and a sample is below `min_fidelity`.
Warns (no fail) when a sample is below `min_fidelity` and not strict.

CLI:
  check_strandedness.py --files <a.strandedness.txt> [<b>...] \
      --declared reverse --min-fidelity 0.70 [--strict]
"""

import argparse
import os
import sys

PE_REV = "1+-,1-+,2++,2--"
PE_FWD = "1++,1--,2+-,2-+"
SE_REV = "+-,-+"
SE_FWD = "++,--"


def _frac(lines, pattern):
    for line in lines:
        if pattern in line:
            return float(line.rsplit(":", 1)[1].strip())
    return None


def parse_one(text):
    """-> dict(sample-agnostic): data_type, reverse, forward, failed."""
    lines = text.splitlines()
    low = text.lower()
    data_type = "paired" if "pairend" in low else ("single" if "singleend" in low else None)
    if data_type == "paired":
        rev, fwd = _frac(lines, PE_REV), _frac(lines, PE_FWD)
    else:
        rev, fwd = _frac(lines, SE_REV), _frac(lines, SE_FWD)
    failed = _frac(lines, "failed to determine")
    return {"data_type": data_type, "reverse": rev, "forward": fwd, "failed": failed}


def evaluate(records, declared, min_fidelity, strict):
    """records: list of (sample, parsed). -> (errors, warnings, summary rows)."""
    errors, warnings, rows = [], [], []
    for sample, p in records:
        if p["data_type"] is None or p["reverse"] is None or p["forward"] is None:
            errors.append("cannot parse orientation for sample %s" % sample)
    if errors:
        return errors, warnings, rows

    dtypes = {p["data_type"] for _, p in records}
    if len(dtypes) > 1:
        errors.append("mixed library types across samples: %s" % ", ".join(sorted(dtypes)))
        return errors, warnings, rows

    for sample, p in records:
        detected = "reverse" if p["reverse"] >= p["forward"] else "forward"
        expected = detected if declared == "unstranded" else declared
        e = p["reverse"] if expected == "reverse" else p["forward"]
        o = p["forward"] if expected == "reverse" else p["reverse"]
        fidelity = e / (e + o) if (e + o) > 0 else None
        rows.append((sample, p["data_type"], detected, declared, e, o, p["failed"], fidelity))

    n_rev = sum(1 for r in rows if r[2] == "reverse")
    n_fwd = len(rows) - n_rev
    majority = "reverse" if n_rev >= n_fwd else "forward"

    if declared == "unstranded":
        errors.append("library_strandedness='unstranded' but data are %s-stranded (majority %s)"
                      % (majority, majority))
        return errors, warnings, rows

    if majority != declared:
        errors.append("declared '%s' but empirical majority is '%s' (%d/%d samples)"
                      % (declared, majority, max(n_rev, n_fwd), len(rows)))

    low = [r[0] for r in rows if r[7] is not None and r[7] < min_fidelity]
    if low:
        msg = "%d sample(s) below fidelity floor %.2f: %s" % (len(low), min_fidelity, ", ".join(low))
        if strict:
            errors.append(msg)
        else:
            warnings.append(msg)
    return errors, warnings, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--files", nargs="+", required=True)
    ap.add_argument("--declared", required=True, choices=["reverse", "forward", "unstranded"])
    ap.add_argument("--min-fidelity", type=float, default=0.70)
    ap.add_argument("--strict", action="store_true")
    ap.add_argument("--summary-out", default="strandedness_summary.tsv")
    args = ap.parse_args()

    records = []
    for f in args.files:
        with open(f) as fh:
            text = fh.read()
        sample = os.path.basename(f).replace(".strandedness.txt", "")
        records.append((sample, parse_one(text)))

    errors, warnings, rows = evaluate(records, args.declared, args.min_fidelity, args.strict)
    for w in warnings:
        print("WARNING CHECK_STRANDEDNESS: %s" % w)
    for e in errors:
        print("ERROR CHECK_STRANDEDNESS: %s" % e, file=sys.stderr)
    if errors:
        sys.exit(1)

    with open(args.summary_out, "w") as out:
        out.write("sample\tdata_type\tdetected_orientation\tdeclared_strandedness"
                  "\texpected_fraction\tcontrary_fraction\tfailed_fraction\tfidelity\n")
        for sample, dtype, detected, declared, e, o, failed, fidelity in sorted(rows):
            out.write("\t".join([sample, dtype, detected, declared,
                                 "%.4f" % e, "%.4f" % o, "%.4f" % failed,
                                 "%.4f" % fidelity]) + "\n")
    fids = [r[7] for r in rows if r[7] is not None]
    median = sorted(fids)[len(fids) // 2] if fids else float("nan")
    print("CHECK_STRANDEDNESS: %s (%d samples, median fidelity %.3f)"
          % (args.declared, len(records), median))


if __name__ == "__main__":
    main()
