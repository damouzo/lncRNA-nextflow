#!/usr/bin/env python3
"""F4.0: validate published deliverables against the machine-readable contract.

Stdlib only (json/csv), no pandas. The contract is `docs/deliverables_contract.json`:
each entry lists required columns and a lower bound on rows. Extra columns are
allowed. Placeholders in the path (`<id>`, `<group>`) are matched by glob, so one
entry covers every contrast/group.

CLI:
  validate_deliverables.py --dir deliverables [--contract docs/deliverables_contract.json] [--describe]
Exits non-zero on any violation.
"""

import argparse
import csv
import glob
import json
import os
import re
import sys

IS_NA = {"", "NA", "NaN", "nan", "None"}


def load_contract(path):
    with open(path) as fh:
        return json.load(fh)


def _read_header(path):
    with open(path) as fh:
        for line in fh:
            if line.strip() == "":
                continue
            return next(csv.reader([line], delimiter="\t"))
    return []


def _count_rows(path):
    n = 0
    with open(path) as fh:
        for i, line in enumerate(fh):
            if line.strip() == "" and i == 0:
                continue
            n += 1
    return max(n - 1, 0)  # minus header


def _type_ok(value, kind):
    if value in IS_NA:
        return True
    if kind in ("string", None):
        return True
    if kind == "bool":
        return value.upper() in ("TRUE", "FALSE", "0", "1")
    if kind == "int":
        try:
            int(value)
            return True
        except ValueError:
            return False
    if kind == "float":
        try:
            float(value)
            return True
        except ValueError:
            return False
    return True


def validate_table(path, spec):
    """-> list of error strings for one table file."""
    errors = []
    if not os.path.exists(path):
        return ["missing file: %s" % path]
    header = _read_header(path)
    required = [c["name"] for c in spec.get("columns", [])]
    missing = [c for c in required if c not in header]
    if missing:
        errors.append("%s: missing columns %s" % (path, ", ".join(missing)))

    min_rows = spec.get("min_rows", 0)
    nrows = _count_rows(path)
    if nrows < min_rows:
        errors.append("%s: %d rows < min_rows %d" % (path, nrows, min_rows))

    idx = {c: header.index(c) for c in required if c in header}
    kinds = {c["name"]: c.get("type") for c in spec.get("columns", [])}
    with open(path) as fh:
        rdr = csv.reader(fh, delimiter="\t")
        next(rdr, None)
        for ln, row in enumerate(rdr, start=2):
            for col, i in idx.items():
                if i >= len(row):
                    errors.append("%s:%d: row shorter than header" % (path, ln))
                    break
                if not _type_ok(row[i], kinds[col]):
                    errors.append("%s:%d: column '%s' value '%s' not %s"
                                  % (path, ln, col, row[i], kinds[col]))
            if len(errors) > 50:
                errors.append("%s: too many errors, truncated" % path)
                return errors
    return errors


def _resolve(base_dir, rel):
    pattern = re.sub(r"<[^>]+>", "*", rel)
    if "*" in pattern:
        return sorted(glob.glob(os.path.join(base_dir, pattern)))
    p = os.path.join(base_dir, rel)
    return [p] if os.path.exists(p) else []


def validate(base_dir, contract):
    errors = []
    for spec in contract.get("files", []):
        matches = _resolve(base_dir, spec["path"])
        if not matches:
            errors.append("no file matches contract path: %s" % spec["path"])
            continue
        kind = spec.get("kind", "table")
        if kind != "table":
            # file/gtf/fasta: existence is the contract.
            for path in matches:
                if not os.path.isfile(path) or os.path.getsize(path) == 0:
                    errors.append("%s: empty or not a file" % path)
            continue
        for path in matches:
            errors.extend(validate_table(path, spec))
    return errors


def describe(contract):
    out = ["Deliverables contract %s v%s" % (contract.get("contract"), contract.get("version"))]
    for spec in contract.get("files", []):
        cols = ", ".join(c["name"] for c in spec.get("columns", []))
        out.append("  %-40s [%s]  %s" % (spec["path"], spec.get("kind", "table"), cols))
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", help="deliverables directory")
    ap.add_argument("--contract", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                       "..", "docs", "deliverables_contract.json"))
    ap.add_argument("--describe", action="store_true")
    args = ap.parse_args()
    contract = load_contract(args.contract)
    if args.describe:
        print(describe(contract))
        return
    if not args.dir:
        ap.error("--dir is required unless --describe is given")
    errors = validate(args.dir, contract)
    if errors:
        for e in errors:
            print("ERROR: %s" % e, file=sys.stderr)
        print("validate_deliverables: FAILED (%d)" % len(errors))
        sys.exit(1)
    print("validate_deliverables: OK")


if __name__ == "__main__":
    main()
