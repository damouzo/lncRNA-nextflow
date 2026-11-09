#!/usr/bin/env python3
"""F2.3: guard the novelty claim.

Without external catalogs there is no evidence, so `novelty_status` is
`not_evaluated` (no failure). With catalogs configured, the guard fails fast when
the evidence is empty or unusable, and otherwise reports an informative verdict
over the novel (MSTRG) *candidates*: the counts below are candidate-scoped (they
include the reference-derived classes c/k/m/n and the ambiguous o that
ANNOTATION_FREEZE later drops); the frozen-catalog novelty figures live in
`deliverables/catalog/catalog_genes.tsv` / `novel_lncrnas.tsv`
(`novelty_class`, `known_catalog_hit`), not here:

  * a BED with 0 features (misparsed / truncated download)
  * a BED that shares no chromosome with the consensus annotation (a genuine
    assembly mismatch; chr1/1 and chrM/MT aliases are normalized first)
  * a catalog run that produced zero known_catalog_hit (novelty would be an
    artefact of the overlap, not a result)
  * a catalog run where the annotated lncRNAs have a 0 hit rate (positive control)

`novelty_status` values:
  * `not_evaluated`          — no catalogs configured
  * `novel_for_all_catalogs` — no novel candidate hits any configured catalog
  * `known`                  — at least one novel candidate hits a catalog

CLI:
  check_novelty.py --overlap catalog_overlap.tsv --consensus-gtf <length-filtered.gtf> \
      [--bed lncipedia=/path/lncipedia.bed ...] --out novelty_status.tsv
"""

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import f1_common as C  # noqa: E402

IS_NA = {"", "NA", "nan", "NaN", "None"}


def read_tsv(path):
    with open(path) as fh:
        rdr = csv.reader(fh, delimiter="\t")
        header = next(rdr, None)
        if header is None:
            return [], []
        return header, [dict(zip(header, row)) for row in rdr if row]


def gtf_chromosomes(path):
    chroms = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t", 1)
            if len(parts) == 2:
                chroms.add(parts[0])
    return chroms


def bed_features_chromosomes(path):
    n = 0
    chroms = set()
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#") or line.startswith("track") \
                    or line.startswith("browser"):
                continue
            n += 1
            chroms.add(line.split("\t", 1)[0])
    return n, chroms


def _is_hit(row):
    return str(row.get("known_catalog_hit", "")).upper() == "TRUE"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--overlap", required=True, help="catalog_overlap.tsv")
    ap.add_argument("--consensus-gtf", required=True, help="consensus/length-filtered GTF")
    ap.add_argument("--bed", action="append", default=[],
                    help="name=path of an external catalog BED (repeatable)")
    ap.add_argument("--out", default="novelty_status.tsv")
    args = ap.parse_args()

    beds = []
    for item in args.bed:
        name, _, path = item.partition("=")
        beds.append((name, path))

    rows = [("n_catalogs", str(len(beds)))]

    if not beds:
        rows.append(("novelty_status", "not_evaluated"))
        rows.append(("novelty_note", "no external catalogs configured"))
        _write(args.out, rows)
        print("check_novelty: no catalogs -> novelty_status=not_evaluated")
        return

    consensus_chroms = gtf_chromosomes(args.consensus_gtf)
    if not consensus_chroms:
        raise SystemExit("consensus GTF has no chromosomes: %s" % args.consensus_gtf)
    # Compare canonical names so a `chr1`/`chrM` catalog agrees with an
    # Ensembl-style `1`/`MT` reference (the overlap step normalizes the same way).
    consensus_canon = {C.canonical_chrom(c) for c in consensus_chroms}

    for name, path in beds:
        if not path or not os.path.exists(path):
            raise SystemExit("catalog BED '%s' not found: %s" % (name, path))
        n, chroms = bed_features_chromosomes(path)
        if n == 0:
            raise SystemExit("catalog BED '%s' has 0 features: %s" % (name, path))
        if not ({C.canonical_chrom(c) for c in chroms} & consensus_canon):
            raise SystemExit(
                "catalog BED '%s' shares no chromosome with the consensus "
                "annotation (assembly mismatch; chr/M/MT aliases are handled)"
                % name)
        rows.append(("n_features_%s" % name, str(n)))

    header, overlap = read_tsv(args.overlap)
    if "known_catalog_hit" not in header:
        raise SystemExit("%s: no 'known_catalog_hit' column" % args.overlap)
    n_total = len(overlap)
    n_hit = sum(1 for r in overlap if _is_hit(r))
    annotated = [r for r in overlap if not str(r.get("ID", "")).startswith("MSTRG")]
    novel = [r for r in overlap if str(r.get("ID", "")).startswith("MSTRG")]
    n_annot_hit = sum(1 for r in annotated if _is_hit(r))
    n_novel_hit = sum(1 for r in novel if _is_hit(r))
    # Counts are over the DISCOVERY CANDIDATE set from CATALOG_OVERLAP, not over the
    # frozen catalog: they include the reference-derived classes (c/k/m/n) and the
    # ambiguous o that ANNOTATION_FREEZE later excludes. The `candidate_` prefix is
    # deliberate so these are never quoted as the deliverable's novel count; the
    # frozen-catalog figures are in catalog/{catalog_genes,novel_lncrnas}.tsv.
    rows.append(("novelty_scope", "candidates"))
    rows.append(("n_candidate_transcripts", str(n_total)))
    rows.append(("n_candidate_known_catalog_hit", str(n_hit)))
    rows.append(("n_candidate_annotated", str(len(annotated))))
    rows.append(("n_candidate_annotated_known_catalog_hit", str(n_annot_hit)))
    rows.append(("n_candidate_novel", str(len(novel))))
    rows.append(("n_candidate_novel_known_catalog_hit", str(n_novel_hit)))

    if n_hit == 0:
        raise SystemExit(
            "catalogs configured but no known_catalog_hit: the novelty claim would "
            "be an artefact of the overlap, not a result")
    if annotated and n_annot_hit == 0:
        raise SystemExit(
            "catalogs configured but annotated lncRNAs have a 0 hit rate "
            "(expected > 0)")

    # Run-wide verdict over the novel (MSTRG) candidates.
    status = "known" if n_novel_hit > 0 else "novel_for_all_catalogs"
    rows.append(("novelty_status", status))
    _write(args.out, rows)
    print("check_novelty: %s (%d/%d novel hit, %d/%d annotated hit)"
          % (status, n_novel_hit, len(novel), n_annot_hit, len(annotated)))


def _write(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["metric", "value"])
        w.writerows(rows)


if __name__ == "__main__":
    main()
