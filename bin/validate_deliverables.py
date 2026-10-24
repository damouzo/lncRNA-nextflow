#!/usr/bin/env python3
"""F4.0: validate published deliverables against the machine-readable contract.

Stdlib only (json/csv), no pandas. The contract is `docs/deliverables_contract.json`:
each entry lists required columns and a lower bound on rows. Extra columns are
allowed. Placeholders in the path (`<id>`, `<group>`) are matched by glob, so one
entry covers every contrast/group.

CLI:
  validate_deliverables.py --dir deliverables [--contract docs/deliverables_contract.json] [--describe]
  validate_deliverables.py --dir deliverables --expect-contrasts comparisons.csv
  validate_deliverables.py --dir deliverables --expect-contrasts contrasts.tsv --contrast-coverage-only
  validate_deliverables.py --expected-ids id1,id2 --found-contrasts id1,id2 --contrast-coverage-only
Exits non-zero on any violation.

`--expect-contrasts` cross-checks the number of `contrasts/<id>/` directories
against the declared contrasts (data rows, header ignored). If the file lists
full `<group>__<contrast>` ids, each one is required to have a directory.
`--contrast-coverage-only` runs just that check, so a final in-run guard can
assert contrast coverage without failing on optional figures or files published
by a later process.

`--expected-ids` / `--found-contrasts` supply the required and delivered
`contrasts/<id>/` ids as comma-separated lists instead of a file and a glob. The
in-run guard uses them because the ids are resolved from the workflow channels:
publishDir copies are asynchronous (dispatched to a thread pool after outputs are
bound), so the published tree cannot be read reliably within the same run.
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
            if not spec.get("optional"):
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


def read_expected_contrasts(path):
    """-> (n_rows, set of full `<group>__<contrast>` ids found in the file).

    Accepts the global contrasts TSV (`name/num/den/batch`, no header in the
    pipeline, header tolerated) or any file listing one contrast per row. A
    declared-name file yields an empty id set (count-only check); a file with
    full ids yields the exact set to require.
    """
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    if not lines:
        return 0, set()
    sep = "\t" if "\t" in lines[0] else ","
    first = [f.strip().lower() for f in lines[0].split(sep)]
    data = lines[1:] if {"contrast_name", "contrast", "numerator"} & set(first) else lines
    ids = set()
    for ln in data:
        # Full ids live in the first (contrast-name) column; derived condition
        # names must never be mistaken for directory ids.
        field = ln.split(sep)[0].strip()
        if "__" in field:
            ids.add(field)
    return len(data), ids


def validate_contrast_coverage(base_dir=None, expected_path=None, found_names=None,
                               expected_ids=None):
    """Compare the declared contrasts against the delivered `contrasts/` folders.

    The declared set comes from `expected_ids` (full `<group>__<contrast>` ids
    resolved by the caller) or, when omitted, from `expected_path` (a file with
    one declared contrast per data row). The delivered set comes from
    `found_names` (ids resolved by the caller) or, when omitted, from globbing
    `<base_dir>/contrasts`. With full ids on both sides each declared id must be
    delivered; with only a declared count, the counts must match.
    """
    errors = []
    if expected_ids is not None:
        expected_ids = set(expected_ids)
        n_expected = len(expected_ids)
    else:
        n_expected, expected_ids = read_expected_contrasts(expected_path)
    if found_names is not None:
        found_ids = set(found_names)
        n_found = len(found_names)
    else:
        cdir = os.path.join(base_dir, "contrasts")
        found = sorted(d for d in glob.glob(os.path.join(cdir, "*"))
                       if os.path.isdir(d)) if os.path.isdir(cdir) else []
        found_ids = {os.path.basename(d) for d in found}
        n_found = len(found)
    if expected_ids:
        missing = sorted(expected_ids - found_ids)
        if missing:
            errors.append("contrast coverage: missing directories for %s"
                          % ", ".join(missing))
    elif n_found != n_expected:
        errors.append("contrast coverage: expected %d, found %d"
                      % (n_expected, n_found))
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
    ap.add_argument("--expect-contrasts",
                    help="file with one data row per declared contrast (header ignored)")
    ap.add_argument("--expected-ids",
                    help="comma-separated <group>__<contrast> ids the run declares; "
                         "each one must be delivered")
    ap.add_argument("--found-contrasts",
                    help="comma-separated ids of the delivered contrasts/ folders; "
                         "used instead of globbing --dir when the caller resolves them")
    ap.add_argument("--contrast-coverage-only", action="store_true",
                    help="only check contrasts/* coverage, skip the full contract")
    ap.add_argument("--describe", action="store_true")
    args = ap.parse_args()
    contract = load_contract(args.contract)
    if args.describe:
        print(describe(contract))
        return
    if not args.dir:
        ap.error("--dir is required unless --describe is given")
    errors = []
    if not args.contrast_coverage_only:
        errors.extend(validate(args.dir, contract))
    if args.expect_contrasts or args.expected_ids:
        found = None
        if args.found_contrasts is not None:
            found = [s for s in args.found_contrasts.split(",") if s]
        expected_ids = None
        if args.expected_ids is not None:
            expected_ids = [s for s in args.expected_ids.split(",") if s]
        errors.extend(validate_contrast_coverage(
            args.dir, args.expect_contrasts, found, expected_ids))
    if errors:
        for e in errors:
            print("ERROR: %s" % e, file=sys.stderr)
        print("validate_deliverables: FAILED (%d)" % len(errors))
        sys.exit(1)
    print("validate_deliverables: OK")


if __name__ == "__main__":
    main()
