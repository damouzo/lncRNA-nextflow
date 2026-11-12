#!/usr/bin/env python3
"""H11 / 0.1: catalog -> DE universe funnel with a reason per dropped gene.

Writes funnel.tsv, a long table (step, gene_id, n, reason) so the run_info tree
shows both the aggregate counts and the per-gene motive for every catalog gene
that never reached the DE universe. The headline case is not_quantified
(transcripts shorter than Salmon's k=31, or contigs missing from genome.fa).

Stdlib only. All inputs except --gene-catalog are optional, so the funnel
degrades gracefully instead of failing when a stage did not run.

usage:
  build_funnel.py --gene-catalog gene_catalog.tsv \\
      --universe MNC_universe.tsv [--universe MS C_universe.tsv ...] \\
      --coding-consensus coding_consensus.tsv --recurrence recurrence_filtered.tsv \\
      --outdir .
"""

import argparse
import csv
import os


def read_tsv(path):
    with open(path) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        return list(rdr)


IS_TRUE = {"TRUE", "true", "1", "yes", "Yes"}


def find_col(header, candidates):
    lowered = {c.lower(): c for c in header}
    for name in candidates:
        if name.lower() in lowered:
            return lowered[name.lower()]
    return None


def unique_genes(rows, col):
    return {r[col] for r in rows if r.get(col)}


# gffcompare classes that are variants of an annotated transcript.
REFERENCE_DERIVED = {"c", "k", "m", "n", "="}


def freeze_class_steps(rec, pcol, ccol, kcol, novel_classes, isoform_classes):
    """Transcripts that pass coding + recurrence, split by what the freeze does.

    frozen_by_class: kept as a catalog entity (class in novel_classes, minus
    isoform_classes). isoform_class: kept in the Salmon index but NOT catalogued
    (class in isoform_classes); its reads are remapped to the reference parent.
    reference_derived_class and excluded_class: passed every filter but left out
    of the freeze; the reason is the class code.
    """
    counts = {}
    for r in rec:
        if str(r.get(pcol, "")) not in IS_TRUE:
            continue
        if kcol and str(r.get(kcol, "")) not in IS_TRUE:
            continue
        code = str(r.get(ccol, "") or "NA")
        counts[code] = counts.get(code, 0) + 1
    out = []
    for code in sorted(counts):
        if code in isoform_classes:
            step = "isoform_class"
        elif code in novel_classes:
            step = "frozen_by_class"
        elif code in REFERENCE_DERIVED:
            step = "reference_derived_class"
        else:
            step = "excluded_class"
        out.append((step, "", counts[code], code))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene-catalog", required=True)
    ap.add_argument("--universe", action="append", default=[],
                    help="<group>_universe.tsv (repeatable; several groups)")
    ap.add_argument("--coding-consensus", default=None)
    ap.add_argument("--recurrence", default=None)
    ap.add_argument("--novel-classes", default=None,
                    help="comma-separated classes ANNOTATION_FREEZE keeps")
    ap.add_argument("--isoform-classes", default=None,
                    help="comma-separated novel-isoform classes kept in the index but not catalogued")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()
    novel_classes = {x for x in (args.novel_classes or "").split(",") if x}
    isoform_classes = {x for x in (args.isoform_classes or "").split(",") if x}

    catalog = read_tsv(args.gene_catalog)
    if catalog:
        id_col = find_col(list(catalog[0].keys()), ["gene_id", "GeneID", "ID", "gene"])
    else:
        id_col = "gene_id"
    if id_col is None:
        raise SystemExit("gene_catalog.tsv has no recognisable gene id column")

    origin_col = find_col(list(catalog[0].keys()), ["origin"]) if catalog else None
    catalog_genes = unique_genes(catalog, id_col)
    novel = {r[id_col] for r in catalog
             if origin_col and str(r.get(origin_col, "")).lower().startswith("novel")}

    universe_genes = set()
    kept_genes = set()
    # Count distinct genes dropped for low expression over the UNION of groups, so
    # the number is on the same scale as de_universe (also a union). Summing the
    # per-group drops inflates it by the number of groups a gene is low in.
    filtered_low_genes = set()
    for path in args.universe:
        rows = read_tsv(path)
        if not rows:
            continue
        gcol = find_col(list(rows[0].keys()), ["gene_id", "GeneID", "ID", "gene"])
        kcol = find_col(list(rows[0].keys()), ["keep"])
        rcol = find_col(list(rows[0].keys()), ["filter_reason", "reason"])
        for r in rows:
            gid = r.get(gcol) if gcol else None
            if not gid:
                continue
            universe_genes.add(gid)
            if kcol and str(r.get(kcol, "")).strip() in IS_TRUE:
                kept_genes.add(gid)
            elif rcol and str(r.get(rcol, "")).strip() == "low_expression":
                filtered_low_genes.add(gid)

    steps = []
    steps.append(("catalog", "", len(catalog_genes), ""))
    steps.append(("catalog_novel", "", len(novel), ""))
    steps.append(("catalog_annotated", "", len(catalog_genes - novel), ""))

    if args.coding_consensus and os.path.exists(args.coding_consensus):
        cons = read_tsv(args.coding_consensus)
        if cons:
            h = list(cons[0].keys())
            cpat = find_col(h, ["CPAT_noncoding"])
            cpc2 = find_col(h, ["CPC2_noncoding"])
            both = find_col(h, ["consensus_noncoding"])
            if cpat:
                steps.append(("cpat_noncoding", "",
                              sum(1 for r in cons if str(r.get(cpat, "")) in IS_TRUE), ""))
            if cpc2:
                steps.append(("cpc2_noncoding", "",
                              sum(1 for r in cons if str(r.get(cpc2, "")) in IS_TRUE), ""))
            if both:
                steps.append(("coding_consensus_noncoding", "",
                              sum(1 for r in cons if str(r.get(both, "")) in IS_TRUE), ""))

    if args.recurrence and os.path.exists(args.recurrence):
        rec = read_tsv(args.recurrence)
        if rec:
            h = list(rec[0].keys())
            pcol = find_col(h, ["passes_recurrence"])
            if pcol:
                steps.append(("recurrence_pass", "",
                              sum(1 for r in rec if str(r.get(pcol, "")) in IS_TRUE), ""))
            ccol = find_col(h, ["class_code"])
            if pcol and ccol and novel_classes:
                steps.extend(freeze_class_steps(rec, pcol, ccol,
                                                find_col(h, ["consensus_noncoding"]),
                                                novel_classes, isoform_classes))

    steps.append(("de_universe", "", len(universe_genes), ""))
    steps.append(("universe_kept", "", len(kept_genes), ""))
    steps.append(("universe_filtered_low_expression", "", len(filtered_low_genes), "low_expression"))
    steps.append(("quantified_catalog", "", len(catalog_genes & universe_genes), ""))

    missing = sorted(catalog_genes - universe_genes)
    steps.append(("not_quantified", "", len(missing), "not_quantified"))
    for gid in missing:
        steps.append(("not_quantified", gid, 1, "not_quantified"))

    os.makedirs(args.outdir, exist_ok=True)
    with open(os.path.join(args.outdir, "funnel.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["step", "gene_id", "n", "reason"])
        w.writerows(steps)

    print("build_funnel: %d catalog genes, %d quantified, %d not_quantified"
          % (len(catalog_genes), len(catalog_genes & universe_genes), len(missing)))


if __name__ == "__main__":
    main()
