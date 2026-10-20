// A10: Orthogonal catalog evidence overlap (two levels)
//
// Thin wrapper: the reciprocal-overlap logic lives in bin/catalog_overlap.py
// (unit-tested with synthetic intervals in bin/tests/test_catalog_overlap.py).
//
//   STRICT (only source of known_catalog_hit / novel figures): overlap / exonic
//     bases of the candidate >= external_catalog_min_reciprocal AND overlap /
//     exonic bases of the external feature >= the same threshold, same strand,
//     evaluated PAIR BY PAIR (overlaps of several features are never summed).
//     Bedtools -f -r convention. Never filters.
//   INFORMATIVE (never filters): <src>_pct = best-pair candidate fraction.
//
// The consensus table is annotated in place: known_catalog_hit, catalog_sources
// and <src>_id/<src>_pct per source (LNCipedia, NONCODE, RNAcentral, FANTOM CAT).
// A missing BED is skipped, never an error. BED6 is matched on its span.

process CATALOG_OVERLAP {

    input:
    path consensus_table
    path consensus_gtf

    output:
    path "catalog_overlap.tsv", emit: overlap_table

    script:
    def bed_args = ['lncipedia', 'noncode', 'rnacentral', 'fantomcat']
        .collect { s -> params."${s}_bed" ? "--${s}-bed ${params."${s}_bed"}" : '' }
        .findAll { it }
        .join(' ')
    """
    catalog_overlap.py \\
        --consensus-table ${consensus_table} \\
        --consensus-gtf ${consensus_gtf} \\
        --min-reciprocal ${params.external_catalog_min_reciprocal} \\
        ${bed_args} \\
        --out catalog_overlap.tsv
    """
}
