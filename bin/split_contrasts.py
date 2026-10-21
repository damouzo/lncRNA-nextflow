#!/usr/bin/env python3
"""Slice D: split the global contrasts file into one single-contrast TSV per
(group, contrast) that actually applies to that norm_group.

Why: REPORT_CONTRAST renders one HTML per contrast and derives the contrasts from
its contrasts file (filtered by the group's conditions). Giving each task a 1-row
file makes it render exactly one report, so publishDir can target
`deliverables/contrasts/<grupo>__<contraste>/` from the task inputs — no R change.

A contrast applies to a group when both numerator and denominator conditions are
present in that group's coldata.

Contrast and group names are validated here because they end up both in the
`<grupo>__<contraste>` directory name and interpolated into the shell scripts of
the downstream modules (build_contrast_tables, make_xlsx, report_contrast).
Allowed charset: [A-Za-z0-9_.-].

CLI:
  split_contrasts.py --coldata coldata_<group>.csv --contrasts contrasts.tsv \
      --group <group> --outdir <dir>
"""

import argparse
import csv
import os
import re

HEADER = ["contrast_name", "numerator", "denominator", "batch"]
SAFE_NAME = re.compile(r"^[A-Za-z0-9_.-]+$")


def check_name(kind, value):
    if "__" in value or "/" in value:
        raise SystemExit(
            "%s '%s' must not contain '__' or '/' "
            "(breaks contrasts/<grupo>__<contraste>/)" % (kind, value))
    if not SAFE_NAME.match(value):
        raise SystemExit(
            "%s '%s' must match [A-Za-z0-9_.-]+ (it is interpolated into shell "
            "commands and directory names)" % (kind, value))


def read_conditions(coldata_csv):
    with open(coldata_csv) as fh:
        rdr = csv.DictReader(fh)
        if "condition" not in (rdr.fieldnames or []):
            raise SystemExit("%s: no 'condition' column" % coldata_csv)
        return {row["condition"] for row in rdr}


def read_contrasts(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            while len(parts) < 4:
                parts.append("")
            rows.append(parts[:4])
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--coldata", required=True)
    ap.add_argument("--contrasts", required=True)
    ap.add_argument("--group", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    conditions = read_conditions(args.coldata)
    os.makedirs(args.outdir, exist_ok=True)
    check_name("group", args.group)

    n = 0
    for name, num, den, batch in read_contrasts(args.contrasts):
        check_name("contrast_name", name)
        if num not in conditions or den not in conditions:
            continue
        out = os.path.join(args.outdir, "%s__%s.tsv" % (args.group, name))
        with open(out, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(HEADER)
            w.writerow([name, num, den, batch])
        n += 1

    print("split_contrasts: %s -> %d contrast file(s) for group %s"
          % (args.contrasts, n, args.group))
    if n == 0:
        print("split_contrasts: WARNING no contrast applies to group %s (conditions: %s)"
              % (args.group, ",".join(sorted(conditions))))


if __name__ == "__main__":
    main()
