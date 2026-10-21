#!/usr/bin/env python3
"""F4.6: assemble the run_info/ audit table set and the run QA report.

Stdlib only. Writes (flat, into --outdir):

  params_used.tsv       key/value of the parameters that shape the run
  versions.tsv          process, tool, version, container — from the versions.yml each
                        tool process emits (reflects the real container, not a probe)
  sample_index.tsv      sample, norm_group, condition, batch
  run_manifest.tsv      key/value provenance (date, pipeline, publish_level, outdir)
  qa_report.tsv         headline QC for reviewers
  external_catalogs.tsv one row per configured external catalog (present true/false)

Versions are best-effort: a tool without a working --version probe records
"not_available". The QA report never gates anything; it only summarises.

CLI:
  build_run_info.py --outdir . --samples samples.tsv --strandedness strandedness_summary.tsv \
      --versions-yml a.yml b.yml --thresholds thresholds_used.tsv --gene-catalog gene_catalog.tsv \
      [--sensitivity tier_final_sensitivity.tsv] --pipeline-version 0.1.0 \
      --publish-level core_raw --outdir-param results [--catalog name=path ...]
"""

import argparse
import csv
import glob
import os
import re
from datetime import datetime, timezone

IS_NA = {"", "NA", "NaN", "nan", "None"}


def read_tsv(path):
    with open(path) as fh:
        rdr = csv.reader(fh, delimiter="\t")
        header = next(rdr, None)
        if header is None:
            return [], []
        return header, [dict(zip(header, row)) for row in rdr if row]


def write_tsv(path, cols, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        w.writerows(rows)


def parse_versions_yml(paths):
    """Collect (process, tool, version, container) from emitted versions.yml files.

    Format (written by bin/emit_versions.sh or the R modules):
        "PROCESS":\n    tool: "x"\n    version: "y"\n    container: "z"

    Every task of a process emits an identical block, so the same row arrives many
    times (one file per task). Keep the first occurrence of each row, in file order,
    so the table reads one row per process.
    """
    rows = []
    seen_files = []
    for pattern in paths:
        seen_files.extend(sorted(glob.glob(pattern)))
    for path in seen_files:
        proc, block = None, {}
        with open(path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line.strip():
                    continue
                if not line[0].isspace() and line.endswith(":"):
                    if proc is not None:
                        rows.append(_version_row(proc, block))
                    proc, block = line[:-1].strip().strip('"'), {}
                elif ":" in line:
                    key, val = line.strip().split(":", 1)
                    block[key.strip()] = val.strip().strip('"')
        if proc is not None:
            rows.append(_version_row(proc, block))

    unique, seen = [], set()
    for row in rows:
        if row not in seen:
            seen.add(row)
            unique.append(row)
    return unique


def _version_row(proc, block):
    return (proc, block.get("tool", ""), block.get("version", ""), block.get("container", ""))


def read_novelty(path):
    if not path or not os.path.exists(path):
        return "not_evaluated"
    _, rows = read_tsv(path)
    for r in rows:
        if r.get("metric") == "novelty_status":
            return r.get("value") or "not_evaluated"
    return "not_evaluated"


def fidelity_range(strandedness_path):
    if not strandedness_path or not os.path.exists(strandedness_path):
        return None
    _, rows = read_tsv(strandedness_path)
    vals = []
    for r in rows:
        v = r.get("fidelity")
        if v is not None and v not in IS_NA:
            try:
                vals.append(float(v))
            except ValueError:
                pass
    if not vals:
        return None
    return min(vals), max(vals)


def thresholds_status(thresholds_path):
    if not thresholds_path or not os.path.exists(thresholds_path):
        return "not_evaluated"
    _, rows = read_tsv(thresholds_path)
    statuses = sorted({r.get("status", "") for r in rows if r.get("status")})
    return ",".join(statuses) if statuses else "not_evaluated"


def count_catalog(gene_catalog):
    if not gene_catalog or not os.path.exists(gene_catalog):
        return None, None
    _, rows = read_tsv(gene_catalog)
    n = len(rows)
    n_tier1 = sum(1 for r in rows if str(r.get("tier_static", "")).upper() in ("TRUE", "1"))
    return n, n_tier1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True, help="task dir where flat run_info files are written")
    ap.add_argument("--samples", required=True)
    ap.add_argument("--strandedness", default=None)
    ap.add_argument("--thresholds", default=None)
    ap.add_argument("--gene-catalog", default=None)
    ap.add_argument("--sensitivity", default=None)
    ap.add_argument("--versions-yml", nargs="+", action="append", default=[],
                    help="versions.yml files emitted by the tool processes. Several "
                         "files per flag and repeated flags both work (argparse would "
                         "otherwise keep only the last flag)")
    ap.add_argument("--novelty-status", default=None,
                    help="CHECK_NOVELTY novelty_status.tsv")
    ap.add_argument("--pipeline-version", default="unknown")
    ap.add_argument("--publish-level", default="core_raw")
    ap.add_argument("--outdir-param", default="results")
    ap.add_argument("--species", default="")
    ap.add_argument("--genome-build", default="")
    ap.add_argument("--library-strandedness", default="")
    ap.add_argument("--design-formula", default="")
    ap.add_argument("--tier1-classes", default="")
    ap.add_argument("--min-transcript-length", default="")
    ap.add_argument("--fdr", default="")
    ap.add_argument("--lfc", default="")
    ap.add_argument("--catalog", action="append", default=[],
                    help="name=path of a configured external catalog (repeatable)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    def out(name):
        return os.path.join(args.outdir, name)

    # params_used.tsv
    params = [
        ("species", args.species), ("genome_build", args.genome_build),
        ("library_strandedness", args.library_strandedness),
        ("design_formula", args.design_formula),
        ("tier1_classes", args.tier1_classes),
        ("min_transcript_length", args.min_transcript_length),
        ("fdr_threshold", args.fdr), ("lfc_threshold", args.lfc),
        ("publish_level", args.publish_level),
    ]
    write_tsv(out("params_used.tsv"), ["param", "value"], params)

    # versions.tsv — one row per emitting process, from the real container
    vrows = parse_versions_yml([p for group in args.versions_yml for p in group])
    if not vrows:
        vrows = [("NA", "NA", "not_available", "none")]
    write_tsv(out("versions.tsv"), ["process", "tool", "version", "container"], vrows)

    # sample_index.tsv
    header, rows = read_tsv(args.samples)
    srows = [[r.get("sample", ""), r.get("norm_group", ""), r.get("condition", ""),
              r.get("batch", "")] for r in rows]
    write_tsv(out("sample_index.tsv"),
              ["sample", "norm_group", "condition", "batch"], srows)

    # run_manifest.tsv
    now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    write_tsv(out("run_manifest.tsv"), ["key", "value"], [
        ("run_date_utc", now),
        ("pipeline_version", args.pipeline_version),
        ("publish_level", args.publish_level),
        ("outdir", args.outdir_param),
        ("n_samples", str(len(rows))),
    ])

    # external_catalogs.tsv
    cat_rows = []
    for item in args.catalog:
        if "=" in item:
            name, path = item.split("=", 1)
        else:
            name, path = item, ""
        present = "TRUE" if path and os.path.exists(path) else "FALSE"
        cat_rows.append([name, path or "NA", present])
    write_tsv(out("external_catalogs.tsv"), ["name", "file", "present"], cat_rows)

    # qa_report.tsv
    n_genes, n_tier1_static = count_catalog(args.gene_catalog)
    fid = fidelity_range(args.strandedness)
    n_tier1_final = ""
    if args.sensitivity and os.path.exists(args.sensitivity):
        _, srows2 = read_tsv(args.sensitivity)
        for r in srows2:
            if r.get("changed") == "baseline":
                n_tier1_final = r.get("n_tier1_final_genes", "")
                break
    qa = [
        ("n_samples", str(len(rows))),
        ("n_catalog_genes", "" if n_genes is None else str(n_genes)),
        ("n_tier1_static_genes", "" if n_tier1_static is None else str(n_tier1_static)),
        ("n_tier1_final_genes", str(n_tier1_final)),
        ("strandedness_fidelity_min", "" if fid is None else "%.4f" % fid[0]),
        ("strandedness_fidelity_max", "" if fid is None else "%.4f" % fid[1]),
        ("thresholds_status", thresholds_status(args.thresholds)),
        ("novelty_status", read_novelty(args.novelty_status)),
    ]
    write_tsv(out("qa_report.tsv"), ["metric", "value"], qa)

    print("build_run_info: wrote params_used/versions/sample_index/run_manifest/"
          "external_catalogs/qa_report into %s" % args.outdir)


if __name__ == "__main__":
    main()
