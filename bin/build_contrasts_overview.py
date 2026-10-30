#!/usr/bin/env python3
"""F4.2/F4.3: one-row-per-contrast overview of the DE deliverables.

Concatenates every per-contrast `summary.tsv` and adds a count of novel Tier-1
genes that were significant in that contrast (join of `de_full_significant.tsv` with
`catalog_genes.tsv`). Stdlib only.

The input files are staged flat with the contrast id encoded in the name
(`<id>__summary.tsv`, `<id>__significant.tsv`) because several share the same
basename and cannot be staged side by side otherwise.

CLI:
  build_contrasts_overview.py --summaries-dir . --catalog catalog_genes.tsv \
      --out contrasts_overview.tsv
"""

import argparse
import csv
import glob
import os

COLS = ["contrast", "n_universe", "n_tested", "n_significant",
        "n_lncrna_full", "n_lncrna_significant",
        "n_novel_tier1_significant", "report"]


def read_tsv(path):
    with open(path) as fh:
        rdr = csv.reader(fh, delimiter="\t")
        header = next(rdr, None)
        if header is None:
            return [], []
        return header, [dict(zip(header, row)) for row in rdr if row]


def _truthy(v):
    return str(v).strip().upper() in ("TRUE", "1", "T", "YES")


def novel_tier1_ids(catalog_path):
    if not os.path.exists(catalog_path):
        raise SystemExit("catalog not found: %s" % catalog_path)
    _, rows = read_tsv(catalog_path)
    return {r.get("gene_id") for r in rows
            if r.get("origin") == "novel" and _truthy(r.get("tier_final"))}


def id_from_name(path, suffix):
    base = os.path.basename(path)
    if not base.endswith(suffix):
        raise SystemExit("unexpected file name: %s" % base)
    return base[: -len(suffix)]


def contrast_from_id(cid):
    """`<group>__<contrast>` -> `<contrast>`. The group is validated upstream to
    never contain '__', so the first separator is the group/contrast boundary."""
    return cid.split("__", 1)[1] if "__" in cid else cid


def build(summaries_dir, catalog_path):
    novel_t1 = novel_tier1_ids(catalog_path)
    summaries = sorted(glob.glob(os.path.join(summaries_dir, "*__summary.tsv")))
    rows = []
    for spath in summaries:
        cid = id_from_name(spath, "__summary.tsv")
        _, srows = read_tsv(spath)
        if not srows:
            raise SystemExit("empty summary: %s" % spath)
        s = srows[0]
        sig_path = os.path.join(summaries_dir, "%s__significant.tsv" % cid)
        n_novel_t1 = 0
        if os.path.exists(sig_path):
            _, sig = read_tsv(sig_path)
            n_novel_t1 = sum(1 for r in sig if r.get("gene_id") in novel_t1)
        rows.append({
            "contrast": cid,
            "n_universe": s.get("n_universe", ""),
            "n_tested": s.get("n_tested", ""),
            "n_significant": s.get("n_significant", ""),
            "n_lncrna_full": s.get("n_lncrna_full", ""),
            "n_lncrna_significant": s.get("n_lncrna_significant", ""),
            "n_novel_tier1_significant": n_novel_t1,
            # REPORT_CONTRAST writes <contrast>_report.html (no group prefix).
            "report": "contrasts/%s/%s_report.html" % (cid, contrast_from_id(cid)),
        })
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summaries-dir", required=True)
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rows = build(args.summaries_dir, args.catalog)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print("build_contrasts_overview: %d contrast(s)" % len(rows))


if __name__ == "__main__":
    main()
