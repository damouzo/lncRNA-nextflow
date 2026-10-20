#!/usr/bin/env python3
"""F2.1: reciprocal exonic overlap between the novel catalog and external catalogs.

Pure Python (stdlib). The logic lives here (not in R) so it is unit-tested with
synthetic intervals; the Nextflow module is a thin wrapper.

Definition (per candidate transcript vs ONE external feature, same strand):
  overlap_bp        = bp shared by the candidate exon union and the feature union
  frac_candidate    = overlap_bp / candidate_exonic_bp
  frac_external     = overlap_bp / external_exonic_bp
  strict  hit  iff frac_candidate >= min_reciprocal AND frac_external >= min_reciprocal
  permissive   iff overlap_bp > 0 (same strand)

The best pair is taken PER transcript and PER source (max min(frac_candidate,
frac_external), tie-broken by overlap_bp). Overlaps of several external features are
NEVER summed: a candidate that only becomes covered by stacking features is not a
strict hit. Sources without exon blocks (BED6) are matched on their span and flagged
`span_only`.

Formats: BED12 (exons from blocks), BED6 (span), GTF (exons per transcript_id).
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import f1_common as C  # noqa: E402

STRANDLESS = {"", ".", "*", "NA", None}


def _int(x):
    try:
        return int(x)
    except ValueError:
        return int(float(x))


def parse_bed(path):
    """-> list of features: {chrom,start,end,strand,name,exons,span_only}.

    BED6 -> span_only (exons = single interval); BED12 -> exon blocks.
    """
    feats = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            chrom, start, end = f[0], _int(f[1]), _int(f[2])
            strand = f[5] if len(f) >= 6 else "."
            name = f[3] if len(f) >= 4 and f[3] else "%s:%d-%d" % (chrom, start, end)
            # BED is 0-based half-open -> convert to 1-based inclusive.
            if len(f) >= 12:
                sizes = [_int(x) for x in f[10].rstrip(",").split(",") if x]
                starts = [_int(x) for x in f[11].rstrip(",").split(",") if x]
                exons = C.merge_intervals([(start + s + 1, start + s + sz)
                                           for s, sz in zip(starts, sizes)])
                feats.append({"chrom": chrom, "start": start + 1, "end": end, "strand": strand,
                              "name": name, "exons": exons, "span_only": False})
            else:
                feats.append({"chrom": chrom, "start": start + 1, "end": end, "strand": strand,
                              "name": name, "exons": [(start + 1, end)], "span_only": True})
    return feats


def parse_gtf_exons(path):
    """-> list of features keyed by transcript_id, exons merged."""
    by_tx = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            attr = C.parse_attributes(f[8])
            tid = attr.get("transcript_id")
            if not tid:
                continue
            d = by_tx.setdefault(tid, {"chrom": f[0], "strand": f[6], "name": tid,
                                       "exons": [], "span_only": False})
            d["exons"].append((_int(f[3]), _int(f[4])))
    feats = []
    for d in by_tx.values():
        ex = C.merge_intervals(d["exons"])
        feats.append({"chrom": d["chrom"], "strand": d["strand"], "name": d["name"],
                      "exons": ex, "start": ex[0][0], "end": ex[-1][1], "span_only": False})
    return feats


def _union_len(ivals):
    return sum(e - s + 1 for s, e in ivals)


def _overlap_bp(a_ivals, b_ivals):
    total = 0
    for s, e in a_ivals:
        for s2, e2 in b_ivals:
            lo, hi = max(s, s2), min(e, e2)
            if hi >= lo:
                total += hi - lo + 1
    return total


def _same_strand(a, b):
    if a in STRANDLESS or b in STRANDLESS:
        return True
    return a == b


def reciprocal(cand_exons, ext_exons, min_reciprocal=0.5):
    """-> (strict, permissive, overlap_bp, frac_candidate, frac_external)."""
    ov = _overlap_bp(cand_exons, ext_exons)
    cb = _union_len(cand_exons) or 1
    eb = _union_len(ext_exons) or 1
    fc, fe = ov / cb, ov / eb
    strict = ov > 0 and fc >= min_reciprocal and fe >= min_reciprocal
    return strict, ov > 0, ov, fc, fe


def annotate_source(candidates, external, source, min_reciprocal=0.5):
    """candidates: {tid: {chrom, strand, exons}}.

    -> {tid: best-pair dict} with keys strict, permissive, frac_candidate,
    frac_external, overlap_bp, ext_id, span_only. Best = max min(frac_c, frac_e).
    """
    by_chrom = {}
    for f in external:
        by_chrom.setdefault(f["chrom"], []).append(f)

    out = {}
    for tid, cand in candidates.items():
        best = None
        for f in by_chrom.get(cand["chrom"], []):
            if not _same_strand(cand.get("strand"), f.get("strand")):
                continue
            strict, perm, ov, fc, fe = reciprocal(cand["exons"], f["exons"], min_reciprocal)
            if not perm:
                continue
            key = (min(fc, fe), ov)
            if best is None or key > best["_key"]:
                best = {"strict": strict, "permissive": True, "frac_candidate": fc,
                        "frac_external": fe, "overlap_bp": ov, "ext_id": f["name"],
                        "span_only": bool(f.get("span_only")), "_key": key}
        if best is None:
            best = {"strict": False, "permissive": False, "frac_candidate": 0.0,
                    "frac_external": 0.0, "overlap_bp": 0, "ext_id": "",
                    "span_only": False, "_key": (0.0, 0)}
        best.pop("_key", None)
        best["source"] = source
        out[tid] = best
    return out


SRCS = ["lncipedia", "noncode", "rnacentral", "fantomcat"]


def _consensus_id_col(header):
    for c in header:
        if c == "ID" or c.lower() in ("transcript_id", "gene_id"):
            return c
    return header[0]


def exons_by_tx(path):
    """transcript_id -> {chrom, strand, exons}. Falls back to the transcript span
    when the GTF carries transcript features without exon records."""
    gtf = C.load_gtf(path, {"transcript", "exon"})
    ex = C.exon_map(gtf)
    out = {}
    for t in gtf["transcript"]:
        tid = t["transcript_id"]
        exons = C.merge_intervals(ex.get(tid, [])) or [(t["start"], t["end"])]
        out[tid] = {"chrom": t["chrom"], "strand": t.get("strand"), "exons": exons}
    return out


def annotate_consensus(header, rows, cands, source_beds, min_reciprocal):
    """Annotate a consensus table with per-source reciprocal overlap.

    Returns the consensus table with the same columns the R module produced,
    appended in the same order: known_catalog_hit, catalog_sources, then
    <src>_id / <src>_pct for each source. <src>_pct is the best-pair candidate
    fraction (informative; never filters). known_catalog_hit = any strict hit.
    """
    id_col = _consensus_id_col(header)
    per = {}
    for s in SRCS:
        p = source_beds.get(s)
        per[s] = annotate_source(cands, parse_bed(p), s, min_reciprocal) if (p and os.path.exists(p)) else {}

    new_header = header + ["known_catalog_hit", "catalog_sources"]
    for s in SRCS:
        new_header += ["%s_id" % s, "%s_pct" % s]

    out = []
    for r in rows:
        tid = r.get(id_col, "")
        cats, vals = [], {}
        for s in SRCS:
            b = per[s].get(tid)
            if b and b["strict"]:
                cats.append(s)
                vals["%s_id" % s] = b["ext_id"]
            else:
                vals["%s_id" % s] = "NA"
            vals["%s_pct" % s] = "%.4f" % (b["frac_candidate"] if b else 0.0)
        r2 = dict(r)
        r2["known_catalog_hit"] = "TRUE" if cats else "FALSE"
        r2["catalog_sources"] = ",".join(cats) if cats else "none"
        r2.update(vals)
        out.append(r2)
    return new_header, out


def _load_candidates(frozen_gtf):
    gtf = C.load_gtf(frozen_gtf, {"transcript", "exon"})
    exons = C.exon_map(gtf)
    cands = {}
    for t in gtf["transcript"]:
        tid = t["transcript_id"]
        ex = C.merge_intervals(exons.get(tid, []))
        if not ex:
            continue
        cands[tid] = {"chrom": t["chrom"], "strand": t.get("strand"), "exons": ex}
    return cands


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-reciprocal", type=float, default=0.5)
    ap.add_argument("--lncipedia-bed")
    ap.add_argument("--noncode-bed")
    ap.add_argument("--rnacentral-bed")
    ap.add_argument("--fantomcat-bed")
    # mode A (drop-in for CATALOG_OVERLAP): annotate the consensus table
    ap.add_argument("--consensus-table")
    ap.add_argument("--consensus-gtf")
    # mode B (hits table): candidates from the frozen GTF
    ap.add_argument("--frozen-gtf")
    args = ap.parse_args()
    beds = {"lncipedia": args.lncipedia_bed, "noncode": args.noncode_bed,
            "rnacentral": args.rnacentral_bed, "fantomcat": args.fantomcat_bed}

    if args.consensus_table and args.consensus_gtf:
        C.ensure_file(args.consensus_table, "consensus table")
        C.ensure_file(args.consensus_gtf, "consensus GTF")
        header, rows = C.read_tsv(args.consensus_table)
        cands = exons_by_tx(args.consensus_gtf)
        new_header, out = annotate_consensus(header, rows, cands, beds, args.min_reciprocal)
        C.write_tsv([[r[c] for c in new_header] for r in out], args.out, new_header)
        n_hit = sum(1 for r in out if r["known_catalog_hit"] == "TRUE")
        print("catalog_overlap (consensus): %d transcripts, %d known_catalog_hit"
              % (len(out), n_hit))
        return

    if not args.frozen_gtf:
        raise SystemExit("catalog_overlap: provide --consensus-table+--consensus-gtf or --frozen-gtf")
    C.ensure_file(args.frozen_gtf, "frozen GTF")
    candidates = _load_candidates(args.frozen_gtf)
    rows = []
    for src, path in beds.items():
        if not path or not os.path.exists(path):
            print("  %s: not configured — skipped" % src)
            continue
        external = parse_bed(path)
        for tid, r in annotate_source(candidates, external, src, args.min_reciprocal).items():
            rows.append((tid, r["source"], "TRUE" if r["strict"] else "FALSE",
                         "TRUE" if r["permissive"] else "FALSE", r["ext_id"],
                         "%.4g" % r["frac_candidate"], "%.4g" % r["frac_external"],
                         r["overlap_bp"], "TRUE" if r["span_only"] else "FALSE"))
    C.write_tsv(rows, args.out,
                ["transcript_id", "source", "strict_hit", "permissive_hit", "external_id",
                 "frac_candidate", "frac_external", "overlap_bp", "span_only"])
    print("catalog_overlap: %d transcripts, %d hits" % (len(candidates), len(rows)))


if __name__ == "__main__":
    main()
