#!/usr/bin/env python3
"""Rewrite the samplesheet's `seq_batch` column from the FASTQ paths.

The value is derived, never hand-written: each FASTQ basename encodes its
flowcell and lane (`<sample>_<lib>_<FLOWCELL>_L<lane>_<read>.fq.gz`). Because a
sample can be sequenced across several lanes (and even different flowcells), a
sample maps to a *set* of flowcell_lane tokens. The pipeline collapses the lanes
of a sample and keeps `extras[0]` (see validate_input.nf), so `seq_batch` is
written as one canonical value per sample: the sorted, unique tokens joined by
`--separator`, repeated on every row of that sample. This keeps the column
deterministic and independent of row order.

`seq_batch` is annotation only (the DE model uses `batch`); it feeds the PCA
point shape and the audit tables.

usage:
  derive_seq_batch.py --samplesheet samplesheet.csv [--out samplesheet.csv]
                      [--separator ';'] [--dry-run]
"""

import argparse
import csv
import os
import re
import sys

# <sample>_<library>_<FLOWCELL>_L<lane>_<read>.fq.gz ; FLOWCELL is the NovaSeq
# token (e.g. 23HJT3LT4), lane is L1/L2/... Read is 1 or 2.
FASTQ_RE = re.compile(r"_([0-9A-Za-z]+)_L(\d+)_[12]\.(?:fq|fastq)\.gz$")


def sniff_delimiter(path):
    return "\t" if path.endswith((".tsv", ".txt")) else ","


def derive_token(fastq_path):
    """Return '<FLOWCELL>_L<lane>' from a FASTQ path, or None if unrecognized."""
    m = FASTQ_RE.search(os.path.basename(fastq_path))
    if not m:
        return None
    return "%s_L%s" % (m.group(1), m.group(2))


def canonical_by_sample(rows, separator):
    """Map each sample to its sorted, unique set of flowcell_lane tokens."""
    per_sample = {}
    for row in rows:
        sample = row.get("sample", "")
        token = derive_token(row.get("fastq_1", ""))
        per_sample.setdefault(sample, set())
        if token:
            per_sample[sample].add(token)
    return {s: separator.join(sorted(t)) for s, t in per_sample.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--out", default=None,
                    help="output samplesheet (default: overwrite the input)")
    ap.add_argument("--separator", default=";",
                    help="joins the flowcell_lane set for multi-lane samples")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    delim = sniff_delimiter(args.samplesheet)
    with open(args.samplesheet, newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        header = list(reader.fieldnames or [])
        rows = list(reader)

    if not header:
        sys.exit("derive_seq_batch: empty samplesheet")
    if "fastq_1" not in header:
        sys.exit("derive_seq_batch: samplesheet has no fastq_1 column")
    if "sample" not in header:
        sys.exit("derive_seq_batch: samplesheet has no sample column")

    # Insert seq_batch after batch (or at the end) when it is not already there.
    if "seq_batch" not in header:
        if "batch" in header:
            header.insert(header.index("batch") + 1, "seq_batch")
        else:
            header.append("seq_batch")

    canonical = canonical_by_sample(rows, args.separator)

    missing = sorted(s for s, t in canonical.items() if not t)
    if missing:
        sys.stderr.write(
            "derive_seq_batch: no flowcell_lane parsed for %d sample(s): %s\n"
            % (len(missing), ", ".join(missing)))

    for row in rows:
        row["seq_batch"] = canonical.get(row.get("sample", ""), "")

    changed = sum(1 for s in canonical if canonical[s])
    print("derive_seq_batch: %d samples, %d with a flowcell_lane token"
          % (len(canonical), changed))

    if args.dry_run:
        for s in sorted(canonical):
            print("  %-32s %s" % (s, canonical[s]))
        return

    out = args.out or args.samplesheet
    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=header, delimiter=delim,
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print("derive_seq_batch: wrote %s" % out)


if __name__ == "__main__":
    main()
