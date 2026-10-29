#!/usr/bin/env python3
"""Positive-control report: where the shRNA knockdowns land in the DE tables.

DDX41 and DHX34 are the knockdown targets of this design, so their shRNA
contrasts are reported with the control gene's log2FoldChange and padj. This is
diagnostic only: a control that does not reach significance is common in noisy
data and does not by itself invalidate the run, so the report never fails the
run. The statuses (pass/fail/not_applicable) are informational.

A contrast is *applicable* to a control gene when its numerator names that gene's
knockdown (``NB4_DDX41_sh1``, ``MSCline_DHX34_sh2``...). Primary-patient contrasts
(MNC/MSC ``*_DDX41_Patient`` vs ``*_Healthy``) are germline heterozygous, not a
knockdown, so they are reported "not_applicable", never "fail".

Writes qc_controls.tsv (one row per applicable contrast x control gene):
  group  contrast  control_gene  control_gene_id  status  log2FoldChange  padj  reason

The DE tables are located through a manifest (group, contrast, file) because
Nextflow stages a `path` input list flattened into the task dir, so the
`<group>_deseq2_results/` folder name is not available at run time.

usage:
  qc_controls.py --manifest qc_de_manifest.tsv --contrasts contrasts.tsv --outdir .
  qc_controls.py --de <path>_DE.tsv ... --contrasts contrasts.tsv --outdir .
"""

import argparse
import csv
import os
import re
import sys

# Ensembl IDs are the stable key: symbols can change between releases, and the DE
# tables are keyed by gene_id.
CONTROLS = [("DDX41", "ENSG00000183258"), ("DHX34", "ENSG00000134815")]


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_contrasts(path):
    """contrast_name -> numerator. Handles both headered (comparisons.csv shape)
    and headerless (the TSV DESEQ2_DE consumes) inputs."""
    out = {}
    with open(path) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if not row or row[0] == "contrast_name":
                continue
            if len(row) >= 2:
                out[row[0]] = row[1]
    return out


def read_manifest(path):
    """(group, contrast, file) rows. Headered if the first cell is 'group'."""
    rows = []
    with open(path) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if not row or row[0] == "group":
                continue
            if len(row) >= 3:
                rows.append((row[0], row[1], row[2]))
    return rows


def contrast_and_group(de_path):
    """Fallback for --de paths: <group>_deseq2_results/<contrast>_DE.tsv."""
    ct = re.sub(r"_DE\.tsv$", "", os.path.basename(de_path))
    parent = os.path.basename(os.path.dirname(de_path.rstrip("/")))
    group = re.sub(r"_deseq2_results$", "", parent) if parent else ""
    return group, ct


def lfc_column(row):
    # Prefer the unshrunken LFC once DESEQ2_DE starts writing it (Tanda 1.2); the
    # shrunken column is sign-equivalent for the < 0 check used here.
    for name in ("lfc_unshrunken", "log2FoldChange"):
        if name in row and row[name] not in ("", "NA", "NaN"):
            return name
    return "log2FoldChange"


def find_gene(rows, gene_id):
    for r in rows:
        if r.get("gene_id") == gene_id:
            return r
    return None


def applicable(numerator, symbol):
    # NB4_DDX41_sh1 / MSCline_DHX34_sh2 -> knockdown; MNC_DDX41_Patient -> not.
    return bool(re.search(r"%s[_-]?sh\d*" % re.escape(symbol), numerator or ""))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--contrasts", required=True)
    ap.add_argument("--manifest", default=None,
                    help="group<TAB>contrast<TAB>file, one row per DE table")
    ap.add_argument("--de", nargs="+", action="append", default=[],
                    help="explicit DE tables (repeated flags and several paths "
                         "per flag both work); used when --manifest is absent")
    ap.add_argument("--de-dir", default=".",
                    help="directory the manifest file names are relative to")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--padj-max", type=float, default=0.05)
    ap.add_argument("--control", action="append", default=[],
                    help="override/adjust control genes as SYMBOL=ENSEMBL_ID")
    args = ap.parse_args()

    controls = list(CONTROLS)
    for item in args.control:
        if "=" in item:
            sym, gid = item.split("=", 1)
            controls = [c for c in controls if c[0] != sym] + [(sym, gid)]

    if args.manifest:
        targets = [(g, ct, os.path.join(args.de_dir, f))
                   for g, ct, f in read_manifest(args.manifest)]
    else:
        targets = [(g, ct, p) for p in (x for group in args.de for x in group)
                   for g, ct in [contrast_and_group(p)]]

    contrasts = read_contrasts(args.contrasts)
    rows = []
    for group, ct, path in targets:
        numerator = contrasts.get(ct, "")
        de = read_tsv(path)
        for symbol, gene_id in controls:
            base = {"group": group, "contrast": ct, "control_gene": symbol,
                    "control_gene_id": gene_id}
            if not applicable(numerator, symbol):
                rows.append({**base, "status": "not_applicable",
                             "log2FoldChange": "", "padj": "",
                             "reason": "numerator is not a %s knockdown" % symbol})
                continue
            rec = find_gene(de, gene_id)
            if rec is None:
                rows.append({**base, "status": "fail", "log2FoldChange": "",
                             "padj": "", "reason": "control gene not in the DE table"})
                continue
            lcol = lfc_column(rec)
            try:
                lfc = float(rec.get(lcol, "nan"))
                padj = float(rec.get("padj", "nan"))
            except ValueError:
                lfc = padj = float("nan")
            ok = (lfc == lfc and padj == padj and lfc < 0.0 and padj < args.padj_max)
            reason = "" if ok else "control did not go down significantly"
            rows.append({**base, "status": "pass" if ok else "fail",
                         "log2FoldChange": "" if lfc != lfc else lfc,
                         "padj": "" if padj != padj else padj,
                         "reason": reason})

    os.makedirs(args.outdir, exist_ok=True)
    cols = ["group", "contrast", "control_gene", "control_gene_id", "status",
            "log2FoldChange", "padj", "reason"]
    with open(os.path.join(args.outdir, "qc_controls.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    fails = [r for r in rows if r["status"] == "fail"]
    n_pass = sum(1 for r in rows if r["status"] == "pass")
    n_na = sum(1 for r in rows if r["status"] == "not_applicable")
    print("qc_controls: %d pass, %d fail, %d not_applicable" % (n_pass, len(fails), n_na))

    if fails:
        msg = "; ".join("%s/%s %s" % (r["group"], r["contrast"], r["control_gene"])
                        for r in fails)
        print("NOTE: positive control did not go down significantly: " + msg,
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
