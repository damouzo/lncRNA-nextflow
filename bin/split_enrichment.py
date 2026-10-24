#!/usr/bin/env python3
"""Repartition a group-level enrichment table into one file per contrast.

`FUNCTIONAL_ENRICHMENT` writes one `ora`/`gsea` file per norm_group, with a
`contrast` column tagging every row. This splits that file by the column so each
contrast folder carries its own enrichment tables, next to the DE tables.

Not a filter that can lose rows: the rows kept per contrast add up to the input
row count, and the header is always written (even for a contrast with no hits), so
an empty result is a valid 0-row table rather than a headerless 0-byte file.

CLI:
  split_enrichment.py --ora <group>_enrichment_ora.tsv --gsea <group>_enrichment_gsea.tsv \
      --contrast <contrast> --outdir <dir>
"""

import argparse
import csv
import os


def read_tsv(path):
    """-> (header, rows). The header is the single source of truth for the output
    columns, so the clusterProfiler layout is not duplicated here."""
    if not os.path.exists(path):
        raise SystemExit("split_enrichment: missing input: %s" % path)
    with open(path) as fh:
        rdr = csv.reader(fh, delimiter="\t")
        header = next(rdr, None)
        if header is None:
            raise SystemExit("split_enrichment: empty input: %s" % path)
        return header, [dict(zip(header, row)) for row in rdr if row]


def split_table(path, contrast, out_path):
    header, rows = read_tsv(path)
    kept = [r for r in rows if r.get("contrast") == contrast]
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(kept)
    return len(kept)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ora", required=True)
    ap.add_argument("--gsea", required=True)
    ap.add_argument("--contrast", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    n_ora = split_table(args.ora, args.contrast,
                        os.path.join(args.outdir, "enrichment_ora.tsv"))
    n_gsea = split_table(args.gsea, args.contrast,
                         os.path.join(args.outdir, "enrichment_gsea.tsv"))
    print("split_enrichment [%s]: %d ORA rows, %d GSEA rows"
          % (args.contrast, n_ora, n_gsea))


if __name__ == "__main__":
    main()
