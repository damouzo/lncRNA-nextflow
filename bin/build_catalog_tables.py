#!/usr/bin/env python3
"""F4.3 / F4.3b: build the deliverables/catalog/ tables.

Stdlib only (csv), no pandas. Inputs are native pipeline artifacts:

  catalog_genes.tsv        gene_catalog.tsv (full universe) LEFT JOIN tier_final_genes.tsv
                           + antisense_status + lnc_position_class (F5.4). One row per gene
                           of the universe; tiers are labels, never a filter.
  catalog_tx.tsv           tier_final_transcripts.tsv (frozen catalog, transcript level)
                           with stable ids resolved from the ASSIGN_CATALOG_IDS GTF
                           (catalog_id, catalog_transcript_id, mstrg_transcript_id and
                           closest_reference_transcript = gffcompare cmp_ref), plus the
                           gene-level annotation (origin, gene_group, lnc_position_class).
  candidates_all.tsv       F4.3b: every merged candidate (class != '=') with
                           filter_status/filter_reason across the discovery funnel
                           (artifact class -> length -> coding -> recurrence -> frozen),
                           plus gene_id and the gene-level annotation so it reads like
                           catalog_genes.tsv.

The candidate funnel replicates TRANSCRIPT_LENGTH_FILTER's logic on the gffcompare
annotated GTF (class_code + summed exon length) and reads the coding/recurrence
decision from recurrence_filtered.tsv (keyed by transcript id). Nothing is dropped.

CLI:
  build_catalog_tables.py --gene-catalog gene_catalog.tsv \
      --tier-final-genes tier_final_genes.tsv --tier-final-transcripts tier_final_transcripts.tsv \
      [--catalog-gtf frozen_lncrna.catalog.gtf] \
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
    "gene_id", "origin", "class_code", "tier", "tier_static", "tier_final",
    "n_tier1_transcripts", "single_exon_low_confidence", "antisense_status",
    "gene_name", "chromosome", "start", "end", "strand", "gene_biotype",
    "parent_biotype", "novel_type", "catalog_id", "locus", "gene_group",
    "ref_gene_id", "ref_gene_name", "ref_gene_relation", "ref_gene_distance",
    "known_catalog_hit", "catalog_sources",
    "x_host_leak_suspect", "thresholds_status", "conservation_mean",
] + LP.POSITION_COLS
TX_COLS = [
    "transcript_id", "gene_id", "class_code", "tier", "tier_static", "tier_final",
    "origin", "gene_group", "lnc_position_class",
    "parent_biotype", "novel_type", "monoexonic", "single_exon_low_confidence",
    "antisense_status", "x_host_leak_suspect", "n_groups_evaluated",
    "thresholds_status",
    "catalog_id", "catalog_transcript_id", "mstrg_transcript_id",
    "closest_reference_transcript",
]
CAND_COLS = [
    "transcript_id", "gene_id", "class_code", "length", "filter_status", "filter_reason",
    "consensus_noncoding", "passes_recurrence", "frozen",
    "catalog_id", "origin", "gene_group", "tier", "lnc_position_class",
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


def gene_group(row):
    """Coarse group from biotype (mirrors build_contrast_tables.gene_group): origin
    alone lumps protein-coding with pseudogenes/TEC/TR under 'reference_other'."""
    origin = row.get("origin", "")
    if origin == "novel":
        return "novel_lncRNA"
    if origin == "annotated_lncrna":
        return "annotated_lncRNA"
    if (row.get("gene_biotype") or "").strip() == "protein_coding":
        return "protein_coding"
    return "other"


def _fill_locus(row):
    """All genes get a locus string; reference rows do not carry one from the
    registry, so it is built from the coordinates already in the catalog."""
    if not _na(row.get("locus")):
        return
    ch, s, e, st = row.get("chromosome"), row.get("start"), row.get("end"), row.get("strand")
    if any(_na(x) for x in (ch, s, e, st)):
        return
    row["locus"] = "%s:%s-%s(%s)" % (ch, s, e, st)


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
        # H26 / 1.7: BUILD_GENE_CATALOG (R) and COMPUTE_TIER_FINAL (Python) compute
        # the same tier rule; assert they agree in-run so a drift aborts the run
        # instead of publishing an inconsistent catalog.
        if g.get("origin") == "novel":
            if not tf:
                raise SystemExit("tier consistency: novel gene %s missing from tier_final_genes.tsv" % gid)
            if "tier_static" in tf and str(g.get("tier_static")) != str(tf.get("tier_static")):
                raise SystemExit("tier consistency: %s tier_static %s (R) != %s (Python)"
                                 % (gid, g.get("tier_static"), tf.get("tier_static")))
            if ("n_tier1_transcripts" in tf and
                    str(g.get("n_tier1_transcripts")) != str(tf.get("n_tier1_transcripts"))):
                raise SystemExit("tier consistency: %s n_tier1_transcripts %s (R) != %s (Python)"
                                 % (gid, g.get("n_tier1_transcripts"), tf.get("n_tier1_transcripts")))
        n_tier1 = g.get("n_tier1_transcripts")
        row = dict(g)
        row["tier_final"] = (tf.get("tier_final") if tf else None) or "FALSE"
        # 1.5: categorical tier {1a,1b,2} for novel genes; NA for reference rows.
        # novel_type already carries the category from BUILD_GENE_CATALOG.
        nt = g.get("novel_type")
        row["tier"] = nt if (g.get("origin") == "novel" and not _na(nt)) else "NA"
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
    # 4b.2b: gene_group + a locus string for every row (the registry only fills
    # the novel ones), so the catalog reads the same way as the DE tables.
    for r in out:
        _fill_locus(r)
        r["gene_group"] = gene_group(r)
    cols = GENES_COLS + [c for c in header if c not in GENES_COLS]
    n = write_tsv(os.path.join(out_dir, "catalog_genes.tsv"), cols, out)
    # F4.3: novel-only view. A filter of the full catalog, not a source of truth:
    # no gene is dropped from catalog_genes.tsv, this is just the novel subset.
    novel = [r for r in out if r.get("origin") == "novel"]
    write_tsv(os.path.join(out_dir, "novel_lncrnas.tsv"), cols, novel)
    n_tier1_genes = sum(1 for r in out if _truthy(r["tier_final"]))
    return n, n_tier1_genes


def parse_catalog_gtf_ids(path):
    """mstrg_transcript_id -> stable ids from the catalog GTF (ASSIGN_CATALOG_IDS).

    The frozen GTF is keyed by the MSTRG (or disambiguated `MSTRG.x__ENST...`)
    transcript id, which is what catalog_tx.tsv carries. The catalog GTF adds the
    LNC gene/transcript ids and keeps gffcompare's `cmp_ref` (closest reference
    transcript), so one pass resolves every stable id without extra joins.
    """
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "transcript":
                continue
            attrs = {}
            for item in parts[8].rstrip(";").split("; "):
                if " " in item:
                    key, val = item.split(" ", 1)
                    attrs[key] = val.strip().strip('"')
            mstrg = attrs.get("mstrg_transcript_id")
            if not mstrg:
                continue
            out[mstrg] = {
                "catalog_id": attrs.get("catalog_id") or "NA",
                "catalog_transcript_id": attrs.get("transcript_id") or "NA",
                "mstrg_transcript_id": mstrg,
                "closest_reference_transcript": attrs.get("cmp_ref") or "NA",
            }
    return out


def load_gene_annotation(path):
    """{gene_id: annotation} from catalog_genes.tsv, for the gen-centric tables.

    Carries the catalog's own column names so `catalog_tx.tsv` and
    `candidates_all.tsv` read like `catalog_genes.tsv` (origin / class_code /
    tier / lnc_position_class / gene_group / catalog_id).
    """
    _, rows = read_tsv(path)
    cols = ("catalog_id", "origin", "class_code", "tier", "lnc_position_class", "gene_group")
    ann = {}
    for r in rows:
        ann[r["gene_id"]] = {c: (r.get(c) or "NA") for c in cols}
    return ann


def build_catalog_tx(tier_final_tx, out_dir, cat_map=None, gene_ann=None):
    header, rows = read_tsv(tier_final_tx)
    for k in ("transcript_id", "gene_id"):
        if k not in header:
            raise SystemExit("%s: missing '%s' column" % (tier_final_tx, k))
    cat_map = cat_map or {}
    gene_ann = gene_ann or {}
    for r in rows:
        # 1.5: transcript-level categorical tier (the frozen catalog is novel only).
        nt = r.get("novel_type")
        r["tier"] = nt if not _na(nt) else "NA"
        # Stable ids from the catalog GTF; all-NA when it is not provided (tests).
        info = cat_map.get(r.get("transcript_id", "")) or {}
        for col in ("catalog_id", "catalog_transcript_id", "mstrg_transcript_id",
                    "closest_reference_transcript"):
            r[col] = info.get(col, "NA")
        # Gene-level annotation joined through gene_id, so the transcript table
        # carries the same classification columns as the gene catalog.
        ann = gene_ann.get(r.get("gene_id", ""), {})
        for col in ("origin", "gene_group", "lnc_position_class"):
            r[col] = ann.get(col, "NA")
        # catalog_id comes from the catalog GTF (transcript-resolved); fall back to
        # the gene-level id when the GTF is absent.
        if r["catalog_id"] == "NA":
            r["catalog_id"] = ann.get("catalog_id", "NA")
    n = write_tsv(os.path.join(out_dir, "catalog_tx.tsv"),
                  TX_COLS + [c for c in header if c not in TX_COLS], rows)
    return n


def parse_gtf_classes_lengths(path):
    """Replicate length_filter.nf: tx -> class_code, tx -> summed exon length,
    tx -> gene_id."""
    tx_class = {}
    tx_len = {}
    tx_gene = {}
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
                gid = attrs.get("gene_id", "").strip('"').strip()
                if gid:
                    tx_gene[tid] = gid
            elif parts[2] == "exon":
                try:
                    tx_len[tid] = tx_len.get(tid, 0) + (int(parts[4]) - int(parts[3]) + 1)
                except ValueError:
                    pass
    return tx_class, tx_len, tx_gene


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
                         length_filtered_gtf=None, gene_ann=None):
    tx_class, tx_len, tx_gene = parse_gtf_classes_lengths(gffcmp_gtf)
    gene_ann = gene_ann or {}
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
        gene = tx_gene.get(tid, "NA")
        ann = gene_ann.get(gene, {})
        rows.append({
            "transcript_id": tid, "gene_id": gene, "class_code": code, "length": length,
            "filter_status": status, "filter_reason": reason,
            "consensus_noncoding": (rec.get(tid, {}) or {}).get("consensus_noncoding", "NA") or "NA",
            "passes_recurrence": (rec.get(tid, {}) or {}).get("passes_recurrence", "NA") or "NA",
            "frozen": "TRUE" if tid in frozen else "FALSE",
            "catalog_id": ann.get("catalog_id", "NA"),
            "origin": ann.get("origin", "NA"),
            "gene_group": ann.get("gene_group", "NA"),
            "tier": ann.get("tier", "NA"),
            "lnc_position_class": ann.get("lnc_position_class", "NA"),
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
    ap.add_argument("--catalog-gtf", default=None,
                    help="ASSIGN_CATALOG_IDS GTF; adds stable LNC ids to catalog_tx.tsv")
    ap.add_argument("--length-filtered-gtf", default=None,
                    help="TRANSCRIPT_LENGTH_FILTER kept GTF; enables parity check")
    ap.add_argument("--recurrence", default=None)
    ap.add_argument("--frozen-gtf", default=None)
    ap.add_argument("--min-length", type=int, default=200)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    tf = index_unique(args.tier_final_genes, "gene_id")
    n_genes, n_tier1 = build_catalog_genes(args.gene_catalog, tf, args.outdir)
    # Gene-level annotation from the catalog just written, so catalog_tx.tsv and
    # candidates_all.tsv carry the same classification columns as catalog_genes.tsv.
    gene_ann = load_gene_annotation(os.path.join(args.outdir, "catalog_genes.tsv"))
    cat_map = parse_catalog_gtf_ids(args.catalog_gtf) if args.catalog_gtf else {}
    n_tx = build_catalog_tx(args.tier_final_transcripts, args.outdir, cat_map, gene_ann)
    print("catalog_genes.tsv: %d genes (%d tier_final)" % (n_genes, n_tier1))
    print("catalog_tx.tsv:    %d transcripts" % n_tx)

    if args.gffcmp_gtf and args.recurrence:
        n_cand, rows = build_candidates_all(args.gffcmp_gtf, args.recurrence,
                                            args.frozen_gtf, args.min_length, args.outdir,
                                            args.length_filtered_gtf, gene_ann)
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
