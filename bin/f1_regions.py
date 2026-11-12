#!/usr/bin/env python3
"""F1 region builder (Phase A, post-freeze).

Builds the regions consumed by FLANK_COVERAGE_RATIO and FLAG_ANTISENSE_LEAKAGE,
for both candidates and the four deterministic control sets. Pure interval logic;
no featureCounts, no Salmon, no DE.

Outputs (in outdir):
  regions.tsv     one row per element: element_id, kind, set, chrom, strand, note
  flank.saf       featureCounts SAF, GeneID = "<element_id>|candidate|left|right"
  antisense.saf   featureCounts SAF, GeneID = "<element_id>|region"

Candidate flanks: mono-exonic frozen transcripts that are Tier 1 (1a, or 1b when
the reference parent gene biotype is in the whitelist). Flanks are +/-flank_bp,
masked against reference-gene exons; for class "i" the flank is the host intron
(coding host only; a non-coding host falls back to genomic flanks).

Controls (deterministic seed, capped):
  pos_flank      annotated mono-exonic expressed lncRNAs
  neg_flank      intronic windows of highly expressed coding genes, length-matched
  pos_antisense  annotated lncRNA / coding-gene opposite-strand exonic overlaps
  null_antisense exons of expressed coding genes with no opposite-strand annotation

Antisense SAF strands are placed on the HOST (coding) strand for candidate x,
pos_antisense and null_antisense alike, so `-s 2` counts host sense reads and
`-s 1` the antisense element: ratio = element/host is comparable to the null and
a low ratio means the antisense signal does not clear the leakage floor.
"""

import argparse
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import f1_common as C  # noqa: E402


def opposite_strand(s):
    return "-" if s == "+" else ("+" if s == "-" else s)


def host_strand(host_id, cand_strand, gene_loc):
    """Strand of the reference host gene; opposite of the candidate as fallback."""
    loc = gene_loc.get(host_id)
    if loc and loc[3] in ("+", "-"):
        return loc[3]
    return opposite_strand(cand_strand)


def genomic_flanks(start, end, flank_bp, chrom_len):
    """Left/right flanks around [start,end], clipped to [1, chrom_len]."""
    out = []
    left = (max(1, start - flank_bp), start - 1)
    right = (end + 1, min(chrom_len, end + flank_bp))
    if left[0] <= left[1]:
        out.append(("left", left))
    if right[0] <= right[1]:
        out.append(("right", right))
    return out


def mask_pieces(interval, masks, min_len):
    """Unmasked pieces of `interval` if total length >= min_len, else None."""
    pieces = C.subtract_intervals(interval, masks)
    pieces = [p for p in pieces if p[1] - p[0] + 1 >= 1]
    if sum(p[1] - p[0] + 1 for p in pieces) < min_len:
        return None
    return pieces


def host_intron(cand_start, cand_end, host_exons):
    """Internal host intron (gap between merged host exons) containing the candidate."""
    merged = C.merge_intervals(host_exons)
    for i in range(len(merged) - 1):
        gap = (merged[i][1] + 1, merged[i + 1][0] - 1)
        if gap[0] <= cand_start and cand_end <= gap[1]:
            return gap
    return None


def antisense_region(tx_exons, host_exons):
    """Merged intersection of a transcript's exons with a host region's exons."""
    if not tx_exons or not host_exons:
        return None
    a = C.merge_intervals(tx_exons)
    b = C.merge_intervals(host_exons)
    out = []
    for x in a:
        for y in b:
            inter = C.interval_intersection(x, y)
            if inter:
                out.append(inter)
    return C.merge_intervals(out) or None


def rank_coding_genes(ref, rec, expressed):
    """Protein-coding genes ranked by best expressed transcript (deterministic)."""
    gb = ref["gene_biotype"]
    best = {}
    for tid, gid in ref["tx2gene"].items():
        if gb.get(gid) != "protein_coding":
            continue
        t = C.strip_version(tid)
        if t in expressed:
            best[gid] = max(best.get(gid, 0), rec.get(t, 0))
    return [g for g, _ in sorted(best.items(), key=lambda kv: (-kv[1], kv[0]))]


def coding_cs_map(ref):
    """{(chrom,strand): merged coding-gene exons}."""
    raw = {}
    for exs in ref["gene_exons_coding"].values():
        for chrom, s, e, strand in exs:
            raw.setdefault((chrom, strand), []).append((s, e))
    return {k: C.merge_intervals(v) for k, v in raw.items()}


def build(ref, frozen, recurrence, outdir, args):
    os.makedirs(outdir, exist_ok=True)
    tier1a = set(args.tier1a_classes)
    tier1 = set(args.tier1_classes)
    parent_bios = set(args.tier1_parent_biotypes)
    annot_bios = set(args.annotated_lncrna_biotypes)

    gene_biotype = ref["gene_biotype"]
    gene_loc = ref["gene_loc"]
    tx_exons_ref = ref["tx_exons"]
    exon_mask = ref["exon_mask"]
    merged_cs = ref["merged_cs"]
    chrom_len = ref["chrom_len"]
    coding_cs = coding_cs_map(ref)
    mask_starts = {ch: C.sorted_starts(v) for ch, v in exon_mask.items()}
    coding_starts = {k: C.sorted_starts(v) for k, v in coding_cs.items()}
    cs_starts = {k: C.sorted_starts(v) for k, v in merged_cs.items()}
    # load_reference already version-strips gene/transcript keys, so this maps
    # cmp_ref (strip_version) and ref_gene_id consistently.
    tx2gene_norm = ref["tx2gene"]

    def outcrop(chrom, strand, starts_map, ivals_map):
        opp = "+" if strand == "-" else "-"
        key = (chrom, opp)
        return ivals_map.get(key, []), starts_map.get(key, [])

    _, rec_rows = C.read_tsv(recurrence)
    rec = {}
    for r in rec_rows:
        tid = C.strip_version(r.get("ID", ""))
        try:
            n = int(float(r.get("n_reps_expressed", "0") or 0))
        except ValueError:
            n = 0
        rec[tid] = n
    expressed = {t for t, n in rec.items() if n >= args.expr_min_reps}

    regions = []
    flank_saf = []
    antisense_saf = []

    # ---- Candidate flanks: mono-exonic Tier 1 frozen transcripts ----
    frozen_exons = C.exon_map(frozen)
    cand_flank = []
    for t in frozen["transcript"]:
        exons = frozen_exons.get(t["transcript_id"], [])
        if len(exons) != 1:
            continue
        parent_id = C.host_gene_of(t, tx2gene_norm)
        pb = gene_biotype.get(parent_id) if parent_id else None
        if C.tier_of(t.get("class_code"), pb, tier1a, tier1, parent_bios) is None:
            continue
        cand_flank.append((t, exons[0], parent_id))

    for t, (start, end), parent_id in cand_flank:
        chrom = t["chrom"]
        clen = chrom_len.get(chrom, end)
        eid = "cand:%s" % t["transcript_id"]
        flank_saf.append(("%s|candidate" % eid, chrom, start, end, t["strand"]))
        regions.append((eid, "flank", "candidate", chrom, t["strand"], "monoexonic"))
        flanks = []
        if t.get("class_code") == "i" and parent_id:
            host_iv = [(s, e) for _, s, e, _ in ref["gene_exons_coding"].get(parent_id, [])]
            hi = host_intron(start, end, host_iv)
            if hi:
                flanks.append(("left", hi))
        if not flanks:
            flanks = genomic_flanks(start, end, args.flank_bp, clen)
        for name, iv in flanks:
            loc = C.local_overlaps(exon_mask.get(chrom, []), mask_starts.get(chrom, []), iv)
            for p in (mask_pieces(iv, loc, args.min_flank_len) or []):
                flank_saf.append(("%s|%s" % (eid, name), chrom, p[0], p[1], t["strand"]))

    # ---- pos_flank: annotated mono-exonic expressed lncRNAs ----
    pos = [tid for tid, exons in tx_exons_ref.items()
           if len(exons) == 1 and C.strip_version(tid) in expressed]
    for tid in C.sample_cap(pos, args.cap, args.seed):
        gid = ref["tx2gene"][tid]
        chrom, _, _, strand = gene_loc[gid]
        start, end = tx_exons_ref[tid][0]
        clen = chrom_len.get(chrom, end)
        eid = "pos_flank:%s" % tid
        flank_saf.append(("%s|candidate" % eid, chrom, start, end, strand))
        regions.append((eid, "flank", "pos_flank", chrom, strand, "annotated_lncrna"))
        for name, iv in genomic_flanks(start, end, args.flank_bp, clen):
            loc = C.local_overlaps(exon_mask.get(chrom, []), mask_starts.get(chrom, []), iv)
            for p in (mask_pieces(iv, loc, args.min_flank_len) or []):
                flank_saf.append(("%s|%s" % (eid, name), chrom, p[0], p[1], strand))

    # ---- neg_flank: intronic windows of highly expressed coding genes ----
    coding_ranked = rank_coding_genes(ref, rec, expressed)
    cand_lengths = [e - s + 1 for _, (s, e), _ in cand_flank] or [args.flank_bp]
    rng = random.Random(args.seed)
    n_neg = 0
    for gid in coding_ranked:
        if n_neg >= args.cap:
            break
        if gid not in gene_loc:
            continue
        chrom, _, _, strand = gene_loc[gid]
        merged = C.merge_intervals([(s, e) for _, s, e, _ in ref["gene_exons_coding"].get(gid, [])])
        if len(merged) < 2:
            continue
        introns = [(merged[i][1] + 1, merged[i + 1][0] - 1)
                   for i in range(len(merged) - 1)]
        introns = [iv for iv in introns if iv[1] - iv[0] + 1 >= 2 * args.min_flank_len]
        if not introns:
            continue
        for k in rng.sample(range(len(introns)), len(introns)):
            if n_neg >= args.cap:
                break
            iv = introns[k]
            length = min(cand_lengths[rng.randrange(len(cand_lengths))], iv[1] - iv[0] + 1)
            w0 = rng.randrange(iv[0], iv[1] - length + 2)
            win = (w0, w0 + length - 1)
            eid = "neg_flank:%s:%d" % (gid, win[0])
            flank_saf.append(("%s|candidate" % eid, chrom, win[0], win[1], strand))
            regions.append((eid, "flank", "neg_flank", chrom, strand, "intronic_window"))
            for name, ivf in genomic_flanks(win[0], win[1], args.flank_bp, chrom_len.get(chrom, win[1])):
                loc = C.local_overlaps(exon_mask.get(chrom, []), mask_starts.get(chrom, []), ivf)
                for p in (mask_pieces(ivf, loc, args.min_flank_len) or []):
                    flank_saf.append(("%s|%s" % (eid, name), chrom, p[0], p[1], strand))
            n_neg += 1

    # ---- Candidate antisense: class "x", plus class "i" on an antisense host ----
    # class i sits fully inside a reference intron; when the host is on the OPPOSITE
    # strand (i_antisense) the stranded protocol resolves it away from the host, so
    # it is treated as an antisense candidate and gets the same exon-overlap leakage
    # control as x. Note the control fires only where the candidate exons overlap the
    # opposite-strand host EXONS (nested/overlapping genes); a purely intronic i has
    # no such overlap and stays antisense_status=not_evaluable. i_sense (host on the
    # same strand) is deliberately excluded: it carries host pre-mRNA risk, which the
    # exon-overlap control cannot capture either.
    for t in frozen["transcript"]:
        cc = t.get("class_code")
        if cc not in ("x", "i"):
            continue
        host = C.host_gene_of(t, tx2gene_norm)
        if cc == "i" and C.intronic_relation(cc, t["strand"], host, gene_loc) != "i_antisense":
            continue
        host_iv = [(s, e) for _, s, e, _ in ref["gene_exons_coding"].get(host, [])]
        reg = antisense_region(frozen_exons.get(t["transcript_id"], []), host_iv)
        if not reg:
            continue
        # The SAF feature is placed on the HOST strand: -s 2 then counts host
        # (sense) reads and -s 1 counts candidate (antisense) reads, so
        # ratio = antisense/sense = candidate/host is directly comparable to the
        # null (a coding gene measured on its own strand).
        hstrand = host_strand(host, t["strand"], gene_loc)
        eid = "cand_antisense:%s" % t["transcript_id"]
        regions.append((eid, "antisense", "candidate", t["chrom"], hstrand, host or ""))
        for p in reg:
            antisense_saf.append(("%s|region" % eid, t["chrom"], p[0], p[1], hstrand))

    # ---- pos_antisense: annotated lncRNA vs coding-gene opposite-strand exons ----
    pos_anti = []
    for tid, exons in tx_exons_ref.items():
        gid = ref["tx2gene"][tid]
        chrom, _, _, strand = gene_loc[gid]
        opp, opp_st = outcrop(chrom, strand, coding_starts, coding_cs)
        span = (min(s for s, _ in exons), max(e for _, e in exons))
        reg = antisense_region(exons, C.local_overlaps(opp, opp_st, span))
        if reg:
            # SAF on the coding (host) strand, same convention as candidates and
            # the null: -s 2 = host sense, -s 1 = annotated lncRNA.
            pos_anti.append((tid, chrom, opposite_strand(strand), reg))
    pos_anti.sort(key=lambda x: x[0])
    if len(pos_anti) > args.cap:
        rng2 = random.Random(args.seed + 1)
        pos_anti = sorted(rng2.sample(pos_anti, args.cap), key=lambda x: x[0])
    for tid, chrom, hstrand, reg in pos_anti:
        eid = "pos_antisense:%s" % tid
        regions.append((eid, "antisense", "pos_antisense", chrom, hstrand, "annotated_pair"))
        for p in reg:
            antisense_saf.append(("%s|region" % eid, chrom, p[0], p[1], hstrand))

    # ---- null_antisense: exons of expressed coding genes, clean opposite strand ----
    n_null = 0
    for gid in coding_ranked:
        if n_null >= args.cap:
            break
        if gid not in gene_loc:
            continue
        chrom, _, _, strand = gene_loc[gid]
        exons = C.merge_intervals([(s, e) for _, s, e, _ in ref["gene_exons_coding"].get(gid, [])])
        if not exons:
            continue
        opp, opp_st = outcrop(chrom, strand, cs_starts, merged_cs)
        clean = [iv for iv in exons if not C.local_overlaps(opp, opp_st, iv)]
        for iv in clean:
            if n_null >= args.cap:
                break
            eid = "null_antisense:%s:%d" % (gid, iv[0])
            regions.append((eid, "antisense", "null_antisense", chrom, strand, ""))
            antisense_saf.append(("%s|region" % eid, chrom, iv[0], iv[1], strand))
            n_null += 1

    C.write_tsv(regions, os.path.join(outdir, "regions.tsv"),
                ["element_id", "kind", "set", "chrom", "strand", "note"])
    C.write_saf(flank_saf, os.path.join(outdir, "flank.saf"))
    C.write_saf(antisense_saf, os.path.join(outdir, "antisense.saf"))

    counts = {}
    for _, kind, set_name, *_ in regions:
        counts[(kind, set_name)] = counts.get((kind, set_name), 0) + 1
    with open(os.path.join(outdir, "regions_stats.tsv"), "w") as fh:
        fh.write("kind\tset\tn_elements\n")
        for (kind, set_name), n in sorted(counts.items()):
            fh.write("%s\t%s\t%d\n" % (kind, set_name, n))
    print("F1 regions: %d elements (%d flank SAF rows, %d antisense SAF rows)"
          % (len(regions), len(flank_saf), len(antisense_saf)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frozen-gtf", required=True)
    ap.add_argument("--reference-gtf", required=True)
    ap.add_argument("--recurrence", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--flank-bp", type=int, default=2000)
    ap.add_argument("--min-flank-len", type=int, default=200)
    ap.add_argument("--cap", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--expr-min-reps", type=int, default=3)
    ap.add_argument("--tier1-classes", default="u,x,i,j")
    ap.add_argument("--tier1a-classes", default="u,x,i")
    ap.add_argument("--tier1-parent-biotypes", default="lncRNA")
    ap.add_argument("--annotated-lncrna-biotypes", default="lncRNA")
    args = ap.parse_args()
    args.tier1_classes = C.parse_list(args.tier1_classes)
    args.tier1a_classes = C.parse_list(args.tier1a_classes)
    args.tier1_parent_biotypes = C.parse_list(args.tier1_parent_biotypes)
    args.annotated_lncrna_biotypes = C.parse_list(args.annotated_lncrna_biotypes)

    C.ensure_file(args.frozen_gtf, "frozen GTF")
    C.ensure_file(args.reference_gtf, "reference GTF")
    C.ensure_file(args.recurrence, "recurrence table")

    ref = C.load_reference(args.reference_gtf, args.annotated_lncrna_biotypes)
    frozen = C.load_gtf(args.frozen_gtf, {"transcript", "exon"})
    if not ref["gene_biotype"]:
        raise SystemExit("ERROR: no 'gene' features in reference GTF")
    if not frozen["transcript"]:
        raise SystemExit("ERROR: no 'transcript' features in frozen GTF")
    build(ref, frozen, args.recurrence, args.outdir, args)


if __name__ == "__main__":
    main()
