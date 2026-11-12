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
        "n_sig_up", "n_sig_down", "n_lncrna_full", "n_lncrna_significant",
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
        n_sig_up = n_sig_down = 0
        if os.path.exists(sig_path):
            _, sig = read_tsv(sig_path)
            n_novel_t1 = sum(1 for r in sig if r.get("gene_id") in novel_t1)
            for r in sig:
                f = _num(r.get("log2FoldChange"))
                if f is None:
                    continue
                if f > 0:
                    n_sig_up += 1
                elif f < 0:
                    n_sig_down += 1
        rows.append({
            "contrast": cid,
            "n_universe": s.get("n_universe", ""),
            "n_tested": s.get("n_tested", ""),
            "n_significant": s.get("n_significant", ""),
            "n_sig_up": n_sig_up,
            "n_sig_down": n_sig_down,
            "n_lncrna_full": s.get("n_lncrna_full", ""),
            "n_lncrna_significant": s.get("n_lncrna_significant", ""),
            "n_novel_tier1_significant": n_novel_t1,
            # REPORT_CONTRAST writes <contrast>_report.html (no group prefix).
            "report": "contrasts/%s/%s_report.html" % (cid, contrast_from_id(cid)),
        })
    return rows


def _num(v):
    if v is None or str(v) in ("", "NA", "NaN", "nan", "None"):
        return None
    try:
        return float(v)
    except ValueError:
        return None


def catalog_annotation(catalog_path):
    """{gene_id: annotation} for the across-contrasts table (from catalog_genes)."""
    _, rows = read_tsv(catalog_path)
    ann = {}
    for r in rows:
        gid = r.get("gene_id")
        if gid is None:
            continue
        cid = r.get("catalog_id")
        ann[gid] = {
            "gene_name": r.get("gene_name") or "NA",
            "gene_group": r.get("gene_group") or "",
            "origin": r.get("origin") or "",
            "tier": r.get("tier") or "",
            "lnc_position_class": r.get("lnc_position_class") or "",
            "gene_id_stable": cid if (cid and str(cid) not in ("", "NA")) else gid,
        }
    return ann


FIXED_ACROSS = ["gene_id", "gene_id_stable", "gene_name", "gene_group", "origin",
                "tier", "lnc_position_class", "n_contrasts_significant"]


def build_across(lncrna_paths, catalog_path):
    """Wide table: one row per lncRNA, one <cid>.log2FoldChange + <cid>.padj block
    per contrast, plus n_contrasts_significant. Column names use the full
    `<group>__<contrast>` id so two groups with the same contrast name do not clash."""
    ann = catalog_annotation(catalog_path)
    per_ct = {}
    order = []
    for path in sorted(lncrna_paths):
        cid = id_from_name(path, "__lncrna_full.tsv")
        order.append(cid)
        _, rows = read_tsv(path)
        m = {}
        for r in rows:
            m[r.get("gene_id")] = (r.get("log2FoldChange", ""), r.get("padj", ""),
                                   _truthy(r.get("significant_padj")))
        per_ct[cid] = m
    genes = set()
    for m in per_ct.values():
        genes |= set(m.keys())
    cols = list(FIXED_ACROSS)
    for cid in order:
        cols += ["%s.log2FoldChange" % cid, "%s.padj" % cid]
    out_rows = []
    for gid in sorted(genes):
        a = ann.get(gid, {})
        row = {
            "gene_id": gid,
            "gene_id_stable": a.get("gene_id_stable", gid),
            "gene_name": a.get("gene_name", "NA"),
            "gene_group": a.get("gene_group", ""),
            "origin": a.get("origin", ""),
            "tier": a.get("tier", ""),
            "lnc_position_class": a.get("lnc_position_class", ""),
        }
        n_sig = 0
        for cid in order:
            lfc, padj, sig = per_ct[cid].get(gid, ("", "", False))
            row["%s.log2FoldChange" % cid] = lfc
            row["%s.padj" % cid] = padj
            if sig:
                n_sig += 1
        row["n_contrasts_significant"] = str(n_sig)
        out_rows.append(row)
    return cols, out_rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summaries-dir", required=True)
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--across-out", help="optional wide lncRNA-across-contrasts table")
    args = ap.parse_args()

    rows = build(args.summaries_dir, args.catalog)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print("build_contrasts_overview: %d contrast(s)" % len(rows))

    if args.across_out:
        lncrna_paths = sorted(glob.glob(os.path.join(args.summaries_dir, "*__lncrna_full.tsv")))
        cols, across = build_across(lncrna_paths, args.catalog)
        with open(args.across_out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            w.writeheader()
            w.writerows(across)
        print("build_contrasts_overview: lncrna_across_contrasts %d row(s) over %d contrast(s)"
              % (len(across), len(lncrna_paths)))


if __name__ == "__main__":
    main()
