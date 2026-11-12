"""Shared helpers for the F1 (tier / flag) analysis.

Pure logic only: GTF parsing, interval arithmetic, SAF writing and deterministic
sampling. No third-party dependencies (runs in the `discovery` container, which
has python3 + featureCounts). Kept deliberately small and heavily unit-tested.
"""

import gzip
import os
import random
import re

try:
    import bisect
except ImportError:  # pragma: no cover
    bisect = None

_ATTR_RE = re.compile(r'(\S+)\s+"([^"]*)"')


def _open(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_attributes(attr):
    """GTF attribute column -> dict. Handles `key "value";` pairs."""
    return {m.group(1): m.group(2) for m in _ATTR_RE.finditer(attr)}


def strip_version(x):
    return re.sub(r"\.\d+$", "", x) if x else x


def canonical_chrom(chrom):
    """Canonical chromosome key for cross-source matching.

    Strips a `chr` prefix and folds mitochondrial aliases (M, MT, chrM, ...) to
    `MT`, so `chr1`/`1` and `chrM`/`MT` compare equal. Any other name is returned
    unchanged: a genuine mismatch must stay a mismatch (the novelty guard relies
    on that to fail fast on a wrong assembly).
    """
    if not chrom:
        return chrom
    c = chrom[3:] if chrom.startswith("chr") else chrom
    if c.lower() in ("m", "mt"):
        return "MT"
    return c


def build_chrom_map(target_chroms):
    """{canonical form -> reference spelling} for a set of reference chromosomes."""
    return {canonical_chrom(t): t for t in target_chroms}


def match_chrom(chrom, chrom_map):
    """Reference spelling for `chrom`, or `chrom` unchanged when it does not map."""
    return chrom_map.get(canonical_chrom(chrom), chrom)


def normalize_chroms(features, chrom_map, key="chrom"):
    """Rewrite `features[*][key]` to the reference spelling. Returns the count changed."""
    changed = 0
    for f in features:
        c = f.get(key)
        t = match_chrom(c, chrom_map)
        if t != c:
            f[key] = t
            changed += 1
    return changed


def parse_list(s):
    return [x.strip() for x in s.split(",") if x.strip()]


def tier_of(class_code, parent_biotype, tier1a, tier1, parent_biotypes):
    """Transcript-level tier from class code and reference parent biotype.

    Single source of truth shared by the region builder (candidate/control
    selection) and the terminal tier computation, so both cannot drift.
    """
    if class_code in tier1a:
        return "1a"
    if class_code in tier1 and parent_biotype in parent_biotypes:
        return "1b"
    return None


def ref_gene_of(rec):
    """Reference parent gene id, preferring `ref_gene_id` over `cmp_ref_gene`.

    Use only when a symbol-or-id value is acceptable. For a real reference-gene
    join (host resolution) use `host_gene_of` with the reference tx2gene map:
    gffcompare's `cmp_ref_gene` holds the reference gene *name* when one exists
    (e.g. "SSU72"), not the Ensembl id, so a direct join against gene_id misses.
    """
    gid = rec.get("ref_gene_id") or rec.get("cmp_ref_gene") or ""
    return strip_version(gid) or None


def host_gene_of(rec, tx2gene=None):
    """Resolve the reference host gene of a gffcompare-annotated transcript.

    `ref_gene_id` (Ensembl id, emitted for re-keyed classes) first; otherwise map
    the reference transcript in `cmp_ref` through `tx2gene` {tid: gid}. Classes
    i/x carry no `ref_gene_id`, and `cmp_ref_gene` is a gene name/symbol, so the
    transcript mapping is the only reliable route for them.
    """
    gid = strip_version(rec.get("ref_gene_id", ""))
    if gid:
        return gid
    tid = strip_version(rec.get("cmp_ref", ""))
    if tx2gene and tid:
        gid = tx2gene.get(tid)
        if gid:
            return strip_version(gid)
    return None


def intronic_relation(class_code, tx_strand, host_gene_id, gene_loc):
    """Sense/antisense split for a class-i transcript.

    class i means "fully contained in a reference intron", and gffcompare does
    NOT require the host to share the query strand: its own source
    (gffcompare.cpp reclass_XStrand) assigns `i` for the opposite strand too. The
    risk is therefore split by host strand:
      i_sense      host on the same strand      -> host pre-mRNA risk
      i_antisense  host on the opposite strand   -> strand resolves it; treat like x
    Returns None for any other class. An unknown host strand is labelled the
    conservative `i_sense` (never silently upgraded to the confident class).
    """
    if class_code != "i":
        return None
    loc = gene_loc.get(host_gene_id) if host_gene_id else None
    if not loc or loc[3] not in ("+", "-"):
        return "i_sense"
    return "i_sense" if loc[3] == tx_strand else "i_antisense"


def load_gtf(path, wanted):
    """Stream a GTF and return {feature_type: [record, ...]}.

    Only feature types listed in `wanted` are kept; a record is a dict with
    chrom/start/end/strand plus every parsed attribute. start/end are 1-based
    inclusive, exactly as in the GTF. For large reference GTFs use
    `load_reference` instead: this keeps one dict per feature and is memory heavy.
    """
    out = {t: [] for t in wanted}
    with _open(path) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            ftype = parts[2]
            if ftype not in wanted:
                continue
            rec = {
                "chrom": parts[0],
                "start": int(parts[3]),
                "end": int(parts[4]),
                "strand": parts[6],
            }
            rec.update(parse_attributes(parts[8]))
            out[ftype].append(rec)
    return out


def load_reference(path, annot_biotypes, coding_biotype="protein_coding"):
    """Memory-lean streaming reference loader for a full Ensembl GTF.

    Returns a compact dict:
      gene_biotype {gid: biotype}
      gene_loc     {gid: (chrom, start, end, strand)}
      tx2gene      {tid: gid}
      tx_exons     {tid: [(start,end)]}  only for annot_biotypes transcripts
      gene_exons_coding {gid: [(chrom,start,end,strand)]}  only coding genes
      exon_mask    {chrom: merged all-gene exons}
      merged_cs    {(chrom,strand): merged all-gene exons}  opposite-strand check
      chrom_len    {chrom: max end}
    Assumes the GTF is hierarchically sorted (gene line before its exons).
    """
    annot_biotypes = set(annot_biotypes)
    gene_biotype = {}
    gene_loc = {}
    tx2gene = {}
    tx_exons = {}
    gene_exons_coding = {}
    cs_raw = {}
    chrom_len = {}
    with _open(path) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            ftype = p[2]
            chrom, strand = p[0], p[6]
            if ftype == "gene":
                a = parse_attributes(p[8])
                gid = strip_version(a.get("gene_id"))
                if gid:
                    gene_biotype[gid] = a.get("gene_biotype")
                    gene_loc[gid] = (chrom, int(p[3]), int(p[4]), strand)
            elif ftype == "transcript":
                a = parse_attributes(p[8])
                tid, gid = strip_version(a.get("transcript_id")), strip_version(a.get("gene_id"))
                if tid and gid:
                    tx2gene[tid] = gid
            elif ftype == "exon":
                a = parse_attributes(p[8])
                tid, gid = strip_version(a.get("transcript_id")), strip_version(a.get("gene_id"))
                start, end = int(p[3]), int(p[4])
                cs_raw.setdefault((chrom, strand), []).append((start, end))
                if end > chrom_len.get(chrom, 0):
                    chrom_len[chrom] = end
                bio = gene_biotype.get(gid)
                if bio == coding_biotype and gid:
                    gene_exons_coding.setdefault(gid, []).append((chrom, start, end, strand))
                if bio in annot_biotypes and tid:
                    tx_exons.setdefault(tid, []).append((start, end))
    merged_cs = {k: merge_intervals(v) for k, v in cs_raw.items()}
    exon_mask = {}
    for (chrom, _strand), ivals in merged_cs.items():
        exon_mask.setdefault(chrom, []).extend(ivals)
    exon_mask = {c: merge_intervals(v) for c, v in exon_mask.items()}
    return {
        "gene_biotype": gene_biotype,
        "gene_loc": gene_loc,
        "tx2gene": tx2gene,
        "tx_exons": tx_exons,
        "gene_exons_coding": gene_exons_coding,
        "exon_mask": exon_mask,
        "merged_cs": merged_cs,
        "chrom_len": chrom_len,
    }


def merge_intervals(intervals):
    """Merge overlapping/adjacent [start, end] (1-based inclusive) intervals."""
    if not intervals:
        return []
    ivals = sorted(intervals)
    merged = [list(ivals[0])]
    for s, e in ivals[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]


def overlaps(a, b):
    """True if 1-based inclusive intervals a and b overlap (any length)."""
    return a[0] <= b[1] and b[0] <= a[1]


def subtract_intervals(target, masks):
    """Return the pieces of `target` not covered by any interval in `masks`.

    Pieces are 1-based inclusive. Handles overlapping masks via merging first.
    """
    s, e = target
    if s > e:
        return []
    merged = merge_intervals([m for m in masks if overlaps(target, m)])
    pieces = []
    cur = s
    for ms, me in merged:
        if ms > cur:
            pieces.append((cur, ms - 1))
        cur = max(cur, me + 1)
        if cur > e:
            break
    if cur <= e:
        pieces.append((cur, e))
    return pieces


def interval_intersection(a, b):
    s = max(a[0], b[0])
    e = min(a[1], b[1])
    return (s, e) if s <= e else None


def sorted_starts(ivals):
    return [iv[0] for iv in ivals]


def local_overlaps(ivals, starts, target):
    """Intervals in `ivals` (sorted by start) overlapping `target`.

    Uses binary search on `starts`, then walks left while the interval can still
    reach the target. Avoids the O(n) filter that dominated the region build.
    """
    s, e = target
    hi = bisect.bisect_right(starts, e)
    out = []
    i = hi - 1
    while i >= 0 and ivals[i][1] >= s:
        out.append(ivals[i])
        i -= 1
    return out


def write_saf(rows, path):
    """Write a featureCounts SAF file. rows: (gene_id, chrom, start, end, strand)."""
    with open(path, "w") as fh:
        fh.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        for gid, chrom, start, end, strand in rows:
            fh.write("%s\t%s\t%d\t%d\t%s\n" % (gid, chrom, start, end, strand))


def write_tsv(rows, path, header):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            # 1.5: a missing value is always NA, never the literal "None".
            fh.write("\t".join("NA" if x is None else str(x) for x in r) + "\n")


def read_tsv(path):
    """Read a headerless or headered TSV -> (header, list_of_dict).

    The first line is treated as the header.
    """
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip()]
    if not lines:
        return [], []
    header = lines[0].split("\t")
    rows = [dict(zip(header, l.split("\t"))) for l in lines[1:]]
    return header, rows


def sample_cap(items, cap, seed):
    """Deterministic cap. Sorts first so the input order never matters."""
    ordered = sorted(items)
    if len(ordered) <= cap:
        return ordered
    rng = random.Random(seed)
    return sorted(rng.sample(ordered, cap))


def exon_map(gtf, chrom=None):
    """{transcript_id: [(start, end), ...]} from exon records."""
    out = {}
    for ex in gtf.get("exon", []):
        if chrom is not None and ex["chrom"] != chrom:
            continue
        out.setdefault(ex["transcript_id"], []).append((ex["start"], ex["end"]))
    return out


def median(xs):
    vals = sorted(x for x in xs if x is not None)
    n = len(vals)
    if n == 0:
        return None
    if n % 2:
        return float(vals[n // 2])
    return (vals[n // 2 - 1] + vals[n // 2]) / 2.0


def percentile(xs, q):
    """Linear-interpolation percentile (q in 0..100). Returns None if empty."""
    vals = sorted(x for x in xs if x is not None)
    if not vals:
        return None
    if len(vals) == 1:
        return float(vals[0])
    pos = (q / 100.0) * (len(vals) - 1)
    lo = int(pos)
    hi = min(lo + 1, len(vals) - 1)
    frac = pos - lo
    return vals[lo] + (vals[hi] - vals[lo]) * frac


def spearman(xs, ys):
    """Spearman rho with average ranks for ties. Returns None if undefined."""
    pairs = [(x, y) for x, y in zip(xs, ys) if x is not None and y is not None]
    n = len(pairs)
    if n < 3:
        return None
    rx = _average_ranks([p[0] for p in pairs])
    ry = _average_ranks([p[1] for p in pairs])
    mx = sum(rx) / n
    my = sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    denx = sum((a - mx) ** 2 for a in rx) ** 0.5
    deny = sum((b - my) ** 2 for b in ry) ** 0.5
    if denx == 0 or deny == 0:
        return None
    return num / (denx * deny)


def _average_ranks(vals):
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    ranks = [0.0] * len(vals)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and vals[order[j + 1]] == vals[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def ensure_file(path, what):
    if not os.path.exists(path):
        raise SystemExit("ERROR: %s not found: %s" % (what, path))
