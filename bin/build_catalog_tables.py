#!/usr/bin/env python3
"""F4.3 / F4.3b: build the deliverables/catalog/ tables.

Stdlib only (csv), no pandas. Inputs are native pipeline artifacts:

  catalog_genes.tsv        gene_catalog.tsv (full universe) LEFT JOIN tier_final_genes.tsv
                           + antisense_status + lnc_position_class (F5.4). One row per gene
                           of the universe; tiers are labels, never a filter.
  catalog_tx.tsv           tier_final_transcripts.tsv (frozen catalog, transcript level).
  candidates_all.tsv       F4.3b: every merged candidate (class != '=') with
                           filter_status/filter_reason across the discovery funnel
                           (artifact class -> length -> coding -> recurrence -> frozen).

The candidate funnel replicates TRANSCRIPT_LENGTH_FILTER's logic on the gffcompare
annotated GTF (class_code + summed exon length) and reads the coding/recurrence
decision from recurrence_filtered.tsv (keyed by transcript id). Nothing is dropped.

CLI:
  build_catalog_tables.py --gene-catalog gene_catalog.tsv \
      --tier-final-genes tier_final_genes.tsv --tier-final-transcripts tier_final_transcripts.tsv \
      [--gffcmp-gtf gffcmp.annotated.gtf] [--length-filtered-gtf transcripts_length_filtered.gtf] \
      [--recurrence recurrence_filtered.tsv] [--frozen-gtf frozen_lncrna.gtf] \
      --min-length 200 --outdir deliverables/catalog
"""

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lnc_position_class as LP  # noqa: E402

IS_NA = {"", "NA", "NaN", "nan", "None", "NA_character_"}

# gffcompare classes kept for discovery vs artifact-prone (mirrors length_filter.nf)
KEEP_CLASSES = {"u", "i", "x", "j", "o", "c", "k", "m", "n"}
EXCLUDE_CLASSES = {"e", "p", "s", "r"}

GENES_COLS = [
    "gene_id", "origin", "class_code", "tier_static", "tier_final",
    "n_tier1_transcripts", "single_exon_low_confidence", "antisense_status",
    "gene_name", "chromosome", "start", "end", "strand", "gene_biotype",
    "parent_biotype", "novel_type", "catalog_id", "locus",
    "ref_gene_id", "ref_gene_name", "ref_gene_relation", "ref_gene_distance",
    "known_catalog_hit", "catalog_sources",
    "x_host_leak_suspect", "thresholds_status", "conservation_mean",
] + LP.POSITION_COLS
TX_COLS = [
    "transcript_id", "gene_id", "class_code", "tier_static", "tier_final",
    "parent_biotype", "novel_type", "monoexonic", "single_exon_low_confidence",
    "antisense_status", "x_host_leak_suspect", "n_groups_evaluated",
    "thresholds_status",
]
CAND_COLS = [
    "transcript_id", "class_code", "length", "filter_status", "filter_reason",
    "consensus_noncoding", "passes_recurrence", "frozen",
]


def read_tsv(path):
    with open(path) as fh:
        rdr = csv.reader(fh, delimiter="\t")
        header = next(rdr, None)
        if header is None:
            return [], []
        return header, [dict(zip(header, row)) for row in rdr if row]


def write_tsv(path, cols, rows):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    return len(rows)


def index_unique(path, key):
    header, rows = read_tsv(path)
    if key not in header:
        raise SystemExit("%s: no '%s' column (have %s)" % (path, key, ",".join(header)))
    out = {}
    for r in rows:
        k = r.get(key)
        if k in out:
            raise SystemExit("duplicate %s '%s' in %s" % (key, k, path))
        out[k] = r
    return out


def _truthy(v):
    return str(v).strip().upper() in ("TRUE", "1", "T", "YES")


def _na(v):
    return v is None or str(v) in IS_NA


def build_catalog_genes(gene_catalog, tier_final, out_dir):
    header, rows = read_tsv(gene_catalog)
    seen = set()
    out = []
    for g in rows:
        gid = g.get("gene_id")
        if gid in seen:
            raise SystemExit("duplicate gene_id '%s' in %s" % (gid, gene_catalog))
        seen.add(gid)
        tf = tier_final.get(gid, {})
        n_tier1 = g.get("n_tier1_transcripts")
        row = dict(g)
        row["tier_final"] = (tf.get("tier_final") if tf else None) or "FALSE"
        row["n_tier1_transcripts"] = "0" if _na(n_tier1) else n_tier1
        row["tier_static"] = g.get("tier_static") if not _na(g.get("tier_static")) else "FALSE"
        row["single_exon_low_confidence"] = (
            g.get("single_exon_low_confidence")
            if not _na(g.get("single_exon_low_confidence")) else "FALSE")
        row["antisense_status"] = (tf.get("antisense_status") if tf else None) or "not_evaluated"
        row["x_host_leak_suspect"] = (tf.get("x_host_leak_suspect") if tf else None) or "FALSE"
        row["thresholds_status"] = (tf.get("thresholds_status") if tf else "") or ""
        out.append(row)
    # F5.4: lncRNA position class vs reference protein-coding genes (gene span).
    LP.annotate_rows(out)
    cols = GENES_COLS + [c for c in header if c not in GENES_COLS]
    n = write_tsv(os.path.join(out_dir, "catalog_genes.tsv"), cols, out)
    # F4.3: novel-only view. A filter of the full catalog, not a source of truth:
    # no gene is dropped from catalog_genes.tsv, this is just the novel subset.
    novel = [r for r in out if r.get("origin") == "novel"]
    write_tsv(os.path.join(out_dir, "novel_lncrnas.tsv"), cols, novel)
    n_tier1_genes = sum(1 for r in out if _truthy(r["tier_final"]))
    return n, n_tier1_genes


def build_catalog_tx(tier_final_tx, out_dir):
    header, rows = read_tsv(tier_final_tx)
    for k in ("transcript_id", "gene_id"):
        if k not in header:
            raise SystemExit("%s: missing '%s' column" % (tier_final_tx, k))
    n = write_tsv(os.path.join(out_dir, "catalog_tx.tsv"),
                  TX_COLS + [c for c in header if c not in TX_COLS], rows)
    return n


def parse_gtf_classes_lengths(path):
    """Replicate length_filter.nf: tx -> class_code, tx -> summed exon length."""
    tx_class = {}
    tx_len = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            attrs = {}
            for item in parts[8].rstrip(";").split("; "):
                if " " in item:
                    key, val = item.split(" ", 1)
                    attrs[key] = val
            tid = attrs.get("transcript_id", "").strip('"').strip()
            if not tid:
                continue
            if parts[2] == "transcript":
                tx_class[tid] = attrs.get("class_code", "").strip('"').strip()
            elif parts[2] == "exon":
                try:
                    tx_len[tid] = tx_len.get(tid, 0) + (int(parts[4]) - int(parts[3]) + 1)
                except ValueError:
                    pass
    return tx_class, tx_len


def count_gtf_transcripts(path):
    ids = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "transcript":
                continue
            for item in parts[8].rstrip(";").split("; "):
                if item.startswith("transcript_id "):
                    ids.add(item.split(" ", 1)[1].strip('"').strip())
                    break
    return len(ids)


def build_candidates_all(gffcmp_gtf, recurrence, frozen_gtf, min_length, out_dir,
                         length_filtered_gtf=None):
    tx_class, tx_len = parse_gtf_classes_lengths(gffcmp_gtf)
    rec = index_unique(recurrence, "ID") if recurrence else {}
    frozen = set()
    if frozen_gtf:
        with open(frozen_gtf) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9 or parts[2] != "transcript":
                    continue
                for item in parts[8].rstrip(";").split("; "):
                    if item.startswith("transcript_id "):
                        frozen.add(item.split(" ", 1)[1].strip('"').strip())
                        break

    rows = []
    for tid, code in tx_class.items():
        if code not in KEEP_CLASSES and code not in EXCLUDE_CLASSES:
            continue
        length = tx_len.get(tid, 0)
        if code in EXCLUDE_CLASSES:
            status, reason = "filtered", "artifact_class"
        elif length < min_length:
            status, reason = "filtered", "short_transcript"
        elif tid not in rec:
            status, reason = "filtered", "not_evaluated"
        elif not _truthy(rec[tid].get("consensus_noncoding")):
            status, reason = "filtered", "coding"
        elif not _truthy(rec[tid].get("passes_recurrence")):
            status, reason = "filtered", "low_recurrence"
        else:
            status, reason = "frozen", ""
        rows.append({
            "transcript_id": tid, "class_code": code, "length": length,
            "filter_status": status, "filter_reason": reason,
            "consensus_noncoding": (rec.get(tid, {}) or {}).get("consensus_noncoding", "NA") or "NA",
            "passes_recurrence": (rec.get(tid, {}) or {}).get("passes_recurrence", "NA") or "NA",
            "frozen": "TRUE" if tid in frozen else "FALSE",
        })
    n = write_tsv(os.path.join(out_dir, "candidates_all.tsv"), CAND_COLS, rows)

    # Parity with the real length filter: the transcripts this replica keeps
    # (class kept + length >= min) must equal TRANSCRIPT_LENGTH_FILTER's kept set.
    if length_filtered_gtf:
        kept = count_gtf_transcripts(length_filtered_gtf)
        mine = sum(1 for r in rows
                   if r["filter_reason"] not in ("artifact_class", "short_transcript"))
        if kept != mine:
            raise SystemExit(
                "candidates_all replica keeps %d transcripts but the length filter "
                "kept %d (%s): class/length logic drifted from length_filter.nf"
                % (mine, kept, length_filtered_gtf))
        print("candidates_all parity: %d kept by replica = %d by length filter" % (mine, kept))
    return n, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene-catalog", required=True)
    ap.add_argument("--tier-final-genes", required=True)
    ap.add_argument("--tier-final-transcripts", required=True)
    ap.add_argument("--gffcmp-gtf", default=None)
    ap.add_argument("--length-filtered-gtf", default=None,
                    help="TRANSCRIPT_LENGTH_FILTER kept GTF; enables parity check")
    ap.add_argument("--recurrence", default=None)
    ap.add_argument("--frozen-gtf", default=None)
    ap.add_argument("--min-length", type=int, default=200)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    tf = index_unique(args.tier_final_genes, "gene_id")
    n_genes, n_tier1 = build_catalog_genes(args.gene_catalog, tf, args.outdir)
    n_tx = build_catalog_tx(args.tier_final_transcripts, args.outdir)
    print("catalog_genes.tsv: %d genes (%d tier_final)" % (n_genes, n_tier1))
    print("catalog_tx.tsv:    %d transcripts" % n_tx)

    if args.gffcmp_gtf and args.recurrence:
        n_cand, rows = build_candidates_all(args.gffcmp_gtf, args.recurrence,
                                            args.frozen_gtf, args.min_length, args.outdir,
                                            args.length_filtered_gtf)
        by_status = {}
        by_reason = {}
        for r in rows:
            by_status[r["filter_status"]] = by_status.get(r["filter_status"], 0) + 1
            by_reason[r["filter_reason"] or "frozen"] = by_reason.get(r["filter_reason"] or "frozen", 0) + 1
        print("candidates_all.tsv: %d candidates" % n_cand)
        print("  by filter_status: %s" % ", ".join("%s=%d" % kv for kv in sorted(by_status.items())))
        print("  by filter_reason: %s" % ", ".join("%s=%d" % kv for kv in sorted(by_reason.items())))
    else:
        print("candidates_all.tsv: skipped (need --gffcmp-gtf and --recurrence)")

    print("build_catalog_tables: OK")


if __name__ == "__main__":
    main()
