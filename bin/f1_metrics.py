#!/usr/bin/env python3
"""F1 metric computation from per-sample featureCounts output (Phase A).

Each counting task writes one single-column featureCounts table per sample; this
script merges the collected tables into matrices and turns them into per-element,
per-norm_group metrics. No thresholds are applied here (see f1_thresholds.py) so a
threshold rule change can be re-run from -resume without re-counting the BAMs.

Subcommands:
  flank      -> flank_ratio_median, opp_strand_ratio_median, n_samples_used, flank_len_used
  antisense  -> region_len, sense_reads, antisense_reads, ratio_median, rho, n

featureCounts SAF GeneID is "<element_id>|<region>" with region in
{candidate, left, right} for flanks and {region} for antisense.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import f1_common as C  # noqa: E402


def _meta_length(field):
    """Meta-feature length. featureCounts writes one number, already the sum of
    the merged segments; a ';' list is summed defensively."""
    if ";" in field:
        return sum(int(x) for x in field.split(";"))
    return int(field)


def parse_featurecounts(path):
    """Parse one featureCounts table -> (meta, counts, cols).

    When a GeneID appears on several SAF rows featureCounts merges it into one
    meta-feature and writes the merged coordinates joined by ';' (e.g.
    "37546084;37548362", also for Chr/Strand). Only Geneid, Length and the count
    columns are meaningful here - the authoritative coordinates live in
    regions.tsv - so Chr/Start/End/Strand are kept verbatim as text and never
    parsed as numbers. Length is the summed meta-feature length, which is what
    the per-bp densities need. If a GeneID still repeats across rows, its length
    and counts are summed so the feature is counted once.
    """
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip()]
    if len(lines) < 2:
        raise SystemExit("ERROR: featureCounts output looks empty: %s" % path)
    if lines[0].startswith("#"):
        header, data = lines[1].split("\t"), lines[2:]
    else:
        header, data = lines[0].split("\t"), lines[1:]
    cols = header[6:]
    meta, counts = {}, {}
    for line in data:
        p = line.split("\t")
        gid = p[0]
        length = _meta_length(p[5])
        vals = [int(x) for x in p[6:]]
        if gid in counts:
            meta[gid]["Length"] += length
            counts[gid] = [a + b for a, b in zip(counts[gid], vals)]
        else:
            meta[gid] = {"Chr": p[1], "Start": p[2], "End": p[3],
                         "Strand": p[4], "Length": length}
            counts[gid] = vals
    return meta, counts, cols


def merge_featurecounts(paths):
    """Merge per-sample single-column tables into one matrix.

    Columns are sorted by BAM basename so -s 2 and -s 1 end up in the same order
    (their per-sample indices are then interchangeable).
    """
    parsed = [parse_featurecounts(p) for p in paths]
    for _, _, cols in parsed:
        if len(cols) != 1:
            raise SystemExit("ERROR: expected one sample column per file, got %d" % len(cols))
    order = sorted(range(len(parsed)), key=lambda i: os.path.basename(parsed[i][2][0]))
    cols = [os.path.basename(parsed[i][2][0]) for i in order]
    meta = parsed[order[0]][0]
    base = set(meta)
    for i in order[1:]:
        if set(parsed[i][0]) != base:
            raise SystemExit("ERROR: featureCounts matrices have different feature sets")
    counts = {gid: [parsed[i][1][gid][0] for i in order] for gid in base}
    return meta, counts, cols


def read_coldata(path):
    """Headerless TSV: sample, norm_group, bam_basename."""
    out = []
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                raise SystemExit("ERROR: coldata needs 3 columns (sample, norm_group, bam): %s" % line)
            out.append((p[0], p[1], os.path.basename(p[2])))
    if not out:
        raise SystemExit("ERROR: coldata is empty: %s" % path)
    return out


def check_assignment(cols, bams):
    """Every matrix column maps to exactly one sample with exactly one norm_group.

    Returns {sample: norm_group}. Fails with the offending name otherwise.
    """
    bam_to_sample = {bam: sample for sample, _g, bam in bams}
    sample_groups = {}
    for sample, group, _bam in bams:
        sample_groups.setdefault(sample, set()).add(group)
    for sample, groups in sample_groups.items():
        if len(groups) > 1:
            raise SystemExit("ERROR: sample '%s' maps to multiple norm_groups: %s"
                             % (sample, sorted(groups)))
    for col in cols:
        b = os.path.basename(col)
        if b not in bam_to_sample:
            raise SystemExit("ERROR: matrix column BAM '%s' has no norm_group in coldata" % b)
    return {sample: next(iter(groups)) for sample, groups in sample_groups.items()}


def read_regions(path):
    _, rows = C.read_tsv(path)
    return {r["element_id"]: r for r in rows}


def column_index(cols, bams):
    """Map each sample to its featureCounts column index (by bam basename)."""
    base_to_idx = {os.path.basename(c): i for i, c in enumerate(cols)}
    idx, missing = {}, []
    for sample, _group, bam in bams:
        if bam in base_to_idx:
            idx[sample] = base_to_idx[bam]
        else:
            missing.append(bam)
    if missing:
        raise SystemExit("ERROR: BAM(s) from coldata not found in featureCounts columns: "
                         + ", ".join(sorted(set(missing))[:5]))
    return idx


def flank_metrics(meta, s2, s1, cols, bams, regions, min_reads):
    check_assignment(cols, bams)
    idx = column_index(cols, bams)
    sample_group = {s: g for s, g, _ in bams}
    out = []
    for eid, info in regions.items():
        if info["kind"] != "flank":
            continue
        cand = "%s|candidate" % eid
        if cand not in s2:
            continue
        cand_len = meta[cand]["Length"]
        flank_ids = [g for g in s2 if g.startswith(eid + "|") and g != cand]
        samples = sorted(idx)
        per_group = {}
        for s in samples:
            j = idx[s]
            c_count = s2[cand][j]
            if c_count < min_reads or cand_len <= 0:
                continue
            dens = []
            for fid in flank_ids:
                flen = meta[fid]["Length"]
                if flen > 0:
                    dens.append(s2[fid][j] / flen)
            if not dens:
                continue
            flank_dens = sum(dens) / len(dens)
            if flank_dens <= 0:
                continue
            fr = (c_count / cand_len) / flank_dens
            opp = s1[cand][j] if cand in s1 else 0
            osr = (opp / c_count) if c_count > 0 else None
            per_group.setdefault(sample_group[s], []).append((fr, osr))
        for group, vals in per_group.items():
            frs = [v[0] for v in vals]
            osrs = [v[1] for v in vals if v[1] is not None]
            out.append((eid, info["set"], group,
                        _fmt(C.median(frs)), _fmt(C.median(osrs)),
                        len(vals), cand_len))
    return out


def antisense_metrics(meta, s2, s1, cols, bams, regions, min_reads):
    check_assignment(cols, bams)
    idx = column_index(cols, bams)
    sample_group = {s: g for s, g, _ in bams}
    out = []
    for eid, info in regions.items():
        if info["kind"] != "antisense":
            continue
        rid = "%s|region" % eid
        if rid not in s2:
            continue
        region_len = meta[rid]["Length"]
        samples = sorted(idx)
        by_group = {}
        for s in samples:
            j = idx[s]
            sense = s2[rid][j]
            if sense < min_reads:
                continue
            anti = s1[rid][j] if rid in s1 else 0
            by_group.setdefault(sample_group[s], []).append((sense, anti))
        for group, vals in by_group.items():
            sense = [v[0] for v in vals]
            anti = [v[1] for v in vals]
            ratios = [a / s for s, a in vals if s > 0]
            rho = C.spearman(sense, anti)
            out.append((eid, info["set"], group,
                        region_len, _fmt(C.median(sense)), _fmt(C.median(anti)),
                        _fmt(C.median(ratios)), _fmt(rho), len(vals)))
    return out


def _fmt(x):
    return "NA" if x is None else ("%.6g" % x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("subcommand", choices=["flank", "antisense"])
    ap.add_argument("--s2", nargs="+", required=True, help="per-sample featureCounts -s 2 tables")
    ap.add_argument("--s1", nargs="+", required=True, help="per-sample featureCounts -s 1 tables")
    ap.add_argument("--regions", required=True)
    ap.add_argument("--coldata", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-reads", type=int, default=10)
    args = ap.parse_args()

    for p in list(args.s2) + list(args.s1) + [args.regions, args.coldata]:
        C.ensure_file(p, "input")
    meta, s2, cols = merge_featurecounts(args.s2)
    _, s1, cols1 = merge_featurecounts(args.s1)
    if cols != cols1:
        raise SystemExit("ERROR: -s 2 and -s 1 matrices have different sample columns")
    bams = read_coldata(args.coldata)
    regions = read_regions(args.regions)

    # Guard against an upstream channel-cardinality bug: the per-sample counters
    # must run once per sample. Refuse to derive thresholds from a partial cohort
    # instead of silently using a single sample.
    n_samples = len({sample for sample, _g, _b in bams})
    if len(args.s2) != n_samples or len(args.s1) != n_samples:
        raise SystemExit(
            "ERROR: incomplete featureCounts input: %d -s2 and %d -s1 tables for "
            "%d coldata samples. The counter ran on the wrong number of samples "
            "(check the BAM/region channel cardinality in workflows/discovery.nf)."
            % (len(args.s2), len(args.s1), n_samples))

    if args.subcommand == "flank":
        rows = flank_metrics(meta, s2, s1, cols, bams, regions, args.min_reads)
        header = ["element_id", "set", "norm_group", "flank_ratio_median",
                  "opp_strand_ratio_median", "n_samples_used", "flank_len_used"]
    else:
        rows = antisense_metrics(meta, s2, s1, cols, bams, regions, args.min_reads)
        header = ["element_id", "set", "norm_group", "region_len", "sense_reads",
                  "antisense_reads", "ratio_median", "rho", "n"]
    C.write_tsv(rows, args.out, header)
    print("F1 %s metrics: %d rows -> %s" % (args.subcommand, len(rows), args.out))


if __name__ == "__main__":
    main()
