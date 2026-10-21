#!/usr/bin/env python3
"""F5.4: lnc_position_class — where each lncRNA sits relative to reference
protein-coding genes (gene span, not exons).

For every lncRNA of the catalog (`origin` in {novel, annotated_lncrna}) the
class is derived from the outer gene span against the reference protein-coding
genes (`origin == reference_other` and `gene_biotype == protein_coding`), which
already live in `catalog_genes.tsv`:

  antisense    overlap with a coding gene on the opposite strand
  intronic     span contained in a coding gene, same strand
  overlapping  partial overlap with a coding gene, same strand
  intergenic   no overlap with any coding gene (lincRNA)

Priority is antisense > intronic > overlapping > intergenic; on a tie the coding
gene with the larger overlap wins, then the smaller gene_id. All overlapping
coding genes are recorded, so ambiguity stays auditable. Nothing is filtered:
`reference_other` rows carry NA.

Library for bin/build_catalog_tables.py; the CLI annotates a catalog table
standalone (audit / tests). TSV I/O is the shared bin/f1_common helper.
"""

import argparse
import bisect
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import f1_common as C  # noqa: E402

LNCRNA_ORIGINS = ("novel", "annotated_lncrna")
CODING_ORIGIN = "reference_other"
CODING_BIOTYPE = "protein_coding"

_PRIORITY = {"antisense": 0, "intronic": 1, "overlapping": 2}

POSITION_COLS = [
    "lnc_position_class",
    "lnc_position_gene_id",
    "lnc_position_overlap_bp",
    "lnc_position_n_genes",
    "lnc_position_genes",
]

_NA = {"", "NA", "NaN", "nan", "None", "NA_character_"}


def _na(v):
    return v is None or str(v).strip() in _NA


def _int(v):
    try:
        return int(str(v).strip())
    except (TypeError, ValueError):
        return None


def classify_pair(lnc, coding):
    """Class and overlap bp of one lncRNA span against one coding gene span.

    Returns (class, overlap_bp) or (None, 0) when the spans do not overlap.
    """
    s = max(lnc["start"], coding["start"])
    e = min(lnc["end"], coding["end"])
    ov = e - s + 1
    if ov <= 0:
        return None, 0
    if lnc["strand"] != coding["strand"]:
        return "antisense", ov
    if lnc["start"] >= coding["start"] and lnc["end"] <= coding["end"]:
        return "intronic", ov
    return "overlapping", ov


def _chr_index(coding_rows):
    """{chrom: (genes, starts, max_end)} with genes sorted by (start, gene_id).

    `max_end[i]` is the running maximum end over genes[0..i]. It is
    non-decreasing, so the first index whose running max reaches `lnc.start` is
    the leftmost gene that can overlap it.
    """
    idx = {}
    for c in coding_rows:
        idx.setdefault(c["chromosome"], []).append(
            (c["start"], c["end"], c["strand"], c["gene_id"]))
    out = {}
    for chrom, genes in idx.items():
        genes.sort(key=lambda t: (t[0], t[3]))
        starts = [g[0] for g in genes]
        max_end = []
        running = -1
        for g in genes:
            running = max(running, g[1])
            max_end.append(running)
        out[chrom] = (genes, starts, max_end)
    return out


def _overlapping(lnc, chrom_idx):
    entry = chrom_idx.get(lnc["chromosome"])
    if not entry:
        return []
    genes, starts, max_end = entry
    # Bounded on both sides: no gene past the first start > lnc.end can overlap,
    # and no gene before the first running-max end >= lnc.start can reach it.
    hi = bisect.bisect_right(starts, lnc["end"])
    lo = bisect.bisect_left(max_end, lnc["start"])
    out = []
    for i in range(lo, hi):
        gs, ge, gstrand, gid = genes[i]
        if ge < lnc["start"]:
            continue
        out.append((gs, ge, gstrand, gid))
    return out


def classify_lnc(lnc, chrom_idx):
    """Position record for one lncRNA against the indexed coding genes."""
    if not lnc.get("chromosome") or _na(lnc["chromosome"]):
        return None
    hits = _overlapping(lnc, chrom_idx)
    if not hits:
        return {"lnc_position_class": "intergenic", "lnc_position_gene_id": "NA",
                "lnc_position_overlap_bp": "0", "lnc_position_n_genes": "0",
                "lnc_position_genes": ""}
    best = None
    best_key = None
    for gs, ge, gstrand, gid in hits:
        cls, ov = classify_pair(lnc, {"start": gs, "end": ge, "strand": gstrand})
        if cls is None:
            continue
        key = (_PRIORITY[cls], -ov, gid)
        if best_key is None or key < best_key:
            best_key, best = key, (cls, gid, ov)
    if best is None:
        # candidate spans touched the window but none overlapped: intergenic
        return {"lnc_position_class": "intergenic", "lnc_position_gene_id": "NA",
                "lnc_position_overlap_bp": "0", "lnc_position_n_genes": "0",
                "lnc_position_genes": ""}
    cls, gid, ov = best
    names = sorted(g[3] for g in hits)
    return {"lnc_position_class": cls, "lnc_position_gene_id": gid,
            "lnc_position_overlap_bp": str(ov), "lnc_position_n_genes": str(len(names)),
            "lnc_position_genes": ",".join(names)}


def _norm_rows(rows):
    """Normalize coordinates to ints and drop rows without them."""
    coding, lncs = [], []
    for r in rows:
        origin = r.get("origin")
        start, end = _int(r.get("start")), _int(r.get("end"))
        if start is None or end is None or not r.get("chromosome"):
            continue
        rec = {"chromosome": r["chromosome"], "start": start, "end": end,
               "strand": r.get("strand"), "gene_id": r.get("gene_id")}
        if origin == CODING_ORIGIN and r.get("gene_biotype") == CODING_BIOTYPE:
            coding.append(rec)
        elif origin in LNCRNA_ORIGINS:
            lncs.append(rec)
    return coding, lncs


def annotate_rows(rows):
    """Add the lnc_position_* keys to every row (in place) and return the rows.

    lncRNAs get a class; `reference_other` rows (and lncRNAs without usable
    coordinates) get NA.
    """
    coding, lncs = _norm_rows(rows)
    chrom_idx = _chr_index(coding)
    by_id = {}
    for lnc in lncs:
        by_id[lnc["gene_id"]] = classify_lnc(lnc, chrom_idx)
    for r in rows:
        rec = by_id.get(r.get("gene_id"))
        if rec is None:
            for col in POSITION_COLS:
                r[col] = "NA"
        else:
            r.update(rec)
    return rows


def main():
    ap = argparse.ArgumentParser(description="F5.4 lncRNA position class")
    ap.add_argument("--catalog", required=True, help="catalog_genes.tsv")
    ap.add_argument("--out", required=True, help="annotated catalog TSV")
    args = ap.parse_args()

    header, rows = C.read_tsv(args.catalog)
    for col in ("gene_id", "origin", "gene_biotype", "chromosome", "start", "end", "strand"):
        if col not in header:
            raise SystemExit("%s: missing '%s' column" % (args.catalog, col))
    annotate_rows(rows)
    out_header = header + [c for c in POSITION_COLS if c not in header]
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    C.write_tsv([[r.get(c, "") for c in out_header] for r in rows], args.out, out_header)
    counts = {}
    for r in rows:
        counts[r["lnc_position_class"]] = counts.get(r["lnc_position_class"], 0) + 1
    print("lnc_position_class: %s" % ", ".join(
        "%s=%d" % kv for kv in sorted(counts.items())))


if __name__ == "__main__":
    sys.exit(main())
