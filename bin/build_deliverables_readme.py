#!/usr/bin/env python3
"""F4.7: generate the deliverables README and its column dictionary from the contract.

The README is derived, never hand-edited: `docs/deliverables_contract.json` is the
single source of truth for the file tree and the required columns. Stdlib only.

Writes two files (names chosen so the Nextflow publishDir can rename them):
  --deliverables-readme   -> deliverables/README.md
  --run-info-readme       -> deliverables/run_info/README.md

CLI:
  build_deliverables_readme.py --contract docs/deliverables_contract.json --outdir . \
      --deliverables-readme __deliverables_README.md \
      --run-info-readme __run_info_README.md [--qa-report qa_report.tsv] \
      [--publish-level core_raw]
"""

import argparse
import csv
import json
import os


def load_contract(path):
    with open(path) as fh:
        return json.load(fh)


def read_kv(path):
    if not path or not os.path.exists(path):
        return []
    with open(path) as fh:
        rdr = csv.reader(fh, delimiter="\t")
        next(rdr, None)
        return [row for row in rdr if row]


def _fmt(value):
    if value is None or value == "":
        return "NA"
    try:
        f = float(value)
    except (TypeError, ValueError):
        return str(value)
    return ("%g" % f)


def build_deliverables_readme(contract, qa_rows, publish_level, fdr=None, lfc=None):
    lines = []
    lines.append("# lncRNA discovery — deliverables")
    lines.append("")
    lines.append("This README and the column dictionary are generated from "
                 "`docs/deliverables_contract.json`; do not edit by hand.")
    lines.append("")
    lines.append("Publish level: **%s**." % publish_level)
    lines.append("")
    lines.append("## Significance")
    lines.append("")
    lines.append("- `significant` = `padj < %s` **and** `|log2FoldChange| >= %s`."
                 % (_fmt(fdr), _fmt(lfc)))
    lines.append("- The cutoff is applied only to this flag, on the **shrunken** "
                 "`log2FoldChange` from `lfcShrink`; it is not passed to DESeq2 "
                 "`results()`, so `pvalue`/`padj` keep their standard meaning.")
    lines.append("- Defaults live in `nextflow.config`/`nextflow_schema.json` "
                 "(`fdr_threshold`, `lfc_threshold`) and are recorded in "
                 "`run_info/params_used.tsv`.")
    lines.append("")
    lines.append("## Opening the TSVs in Excel")
    lines.append("")
    lines.append("Do **not** double-click the `.tsv` files. Use **Data → From Text/CSV** "
                 "(or *Importar desde texto*) and set the columns holding gene symbols "
                 "to **Text**: Excel silently turns `SEPT2` into a date and `MARCH1` into "
                 "`1-Mar`. The same applies to ids that look numeric.")
    lines.append("")
    lines.append("## Honesty note")
    lines.append("")
    lines.append("- Tiers are **labels, not filters**. No table drops a gene by tier or "
                 "significance; every gene keeps its marks (`single_exon_low_confidence`, "
                 "`antisense_status`, `x_host_leak_suspect`).")
    lines.append("- `tier_final` equals `tier_static` in this version (the antisense flag "
                 "annotates but does not demote).")
    lines.append("- `NA` means *not evaluated*, never *clean*.")
    lines.append("")
    lines.append("## Files")
    lines.append("")
    for spec in contract.get("files", []):
        cols = ", ".join(c["name"] for c in spec.get("columns", []))
        lines.append("### `%s`" % spec["path"])
        lines.append("")
        lines.append("- kind: %s" % spec.get("kind", "table"))
        lines.append("- min rows: %s" % spec.get("min_rows", 0))
        if cols:
            lines.append("- required columns: %s" % cols)
        lines.append("")
    if qa_rows:
        lines.append("## Run summary")
        lines.append("")
        lines.append("| metric | value |")
        lines.append("|---|---|")
        for row in qa_rows:
            if len(row) >= 2:
                lines.append("| %s | %s |" % (row[0], row[1]))
        lines.append("")
    lines.append("## Provenance")
    lines.append("")
    lines.append("Run metadata, thresholds, strandedness and external catalogs live under "
                 "`../run_info/` (audit) and `run_info/` (summary). The execution trace, "
                 "report, timeline and DAG are under `run_info/pipeline_info/`.")
    lines.append("")
    return "\n".join(lines)


def build_run_info_readme(contract, qa_rows, publish_level):
    lines = []
    lines.append("# run_info — summary")
    lines.append("")
    lines.append("Curated metadata for this run (publish level **%s**). Generated from "
                 "`docs/deliverables_contract.json`." % publish_level)
    lines.append("")
    lines.append("| file | content |")
    lines.append("|---|---|")
    lines.append("| `versions.tsv` | pipeline and tool versions (best effort) |")
    lines.append("| `sample_index.tsv` | sample, norm_group, condition, batch |")
    lines.append("| `qa_report.tsv` | headline QC; `NA` = not evaluated |")
    lines.append("| `pipeline_info/` | Nextflow trace, report, timeline, DAG |")
    lines.append("")
    if qa_rows:
        lines.append("| metric | value |")
        lines.append("|---|---|")
        for row in qa_rows:
            if len(row) >= 2:
                lines.append("| %s | %s |" % (row[0], row[1]))
        lines.append("")
    lines.append("The full audit set (parameters, thresholds, strandedness summary, "
                 "external catalogs, manifest) lives at `<outdir>/run_info/`.")
    lines.append("")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--contract", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--deliverables-readme", default="__deliverables_README.md")
    ap.add_argument("--run-info-readme", default="__run_info_README.md")
    ap.add_argument("--qa-report", default=None)
    ap.add_argument("--publish-level", default="core_raw")
    ap.add_argument("--fdr", default=None)
    ap.add_argument("--lfc", default=None)
    args = ap.parse_args()

    contract = load_contract(args.contract)
    qa_rows = read_kv(args.qa_report)
    os.makedirs(args.outdir, exist_ok=True)

    with open(os.path.join(args.outdir, args.deliverables_readme), "w") as fh:
        fh.write(build_deliverables_readme(contract, qa_rows, args.publish_level,
                                           fdr=args.fdr, lfc=args.lfc))
    with open(os.path.join(args.outdir, args.run_info_readme), "w") as fh:
        fh.write(build_run_info_readme(contract, qa_rows, args.publish_level))
    print("build_deliverables_readme: wrote %s and %s"
          % (args.deliverables_readme, args.run_info_readme))


if __name__ == "__main__":
    main()
