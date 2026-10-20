#!/usr/bin/env python3
"""Consistency check between BUILD_GENE_CATALOG and COMPUTE_TIER_FINAL.

Both processes compute the same Tier 1 rule (class tier AND mono-exonic rule) for
novel genes: `BUILD_GENE_CATALOG` in R writes `gene_catalog.tsv`, and
`COMPUTE_TIER_FINAL` in Python writes `tier_final_genes.tsv`. This script asserts
they agree on `tier_static` and `n_tier1_transcripts` for every novel gene, so a
drift in either implementation is caught after a run.

Run after a run:
  python3 bin/tests/check_tier_consistency.py \
      results/discovery/gene_catalog/gene_catalog.tsv \
      results/f1/tier_final/tier_final_genes.tsv

Exits non-zero on any mismatch. No third-party dependencies.
"""

import sys

sys.path.insert(0, __file__.rsplit("/", 2)[0])
import f1_common as C  # noqa: E402


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: check_tier_consistency.py gene_catalog.tsv tier_final_genes.tsv")
    gc_path, tf_path = sys.argv[1], sys.argv[2]
    _, gc_rows = C.read_tsv(gc_path)
    _, tf_rows = C.read_tsv(tf_path)

    gc = {r["gene_id"]: r for r in gc_rows if r.get("origin") == "novel"}
    tf = {r["gene_id"]: r for r in tf_rows}

    mismatches = []
    for gid, g in gc.items():
        t = tf.get(gid)
        if t is None:
            mismatches.append((gid, "missing from tier_final_genes.tsv"))
            continue
        if g.get("tier_static") != t.get("tier_static"):
            mismatches.append((gid, "tier_static %s != %s" % (g.get("tier_static"), t.get("tier_static"))))
        if str(g.get("n_tier1_transcripts")) != str(t.get("n_tier1_transcripts")):
            mismatches.append((gid, "n_tier1_transcripts %s != %s"
                               % (g.get("n_tier1_transcripts"), t.get("n_tier1_transcripts"))))

    print("novel genes: gene_catalog=%d, tier_final=%d, mismatches=%d"
          % (len(gc), len(tf), len(mismatches)))
    for gid, msg in mismatches[:50]:
        print("  MISMATCH %s: %s" % (gid, msg))
    if mismatches:
        sys.exit(1)
    print("tier_static consistent between BUILD_GENE_CATALOG and COMPUTE_TIER_FINAL")


if __name__ == "__main__":
    main()
