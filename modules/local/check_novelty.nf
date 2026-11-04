// F2.3: novelty guard. Cheap, always runs: with no catalogs it writes
// novelty_status=not_evaluated; with catalogs it fails fast on an empty or
// mismatched BED (0 features, or no chromosome in common with the consensus)
// and otherwise reports known / novel_for_all_catalogs over the novel MSTRG set.

process CHECK_NOVELTY {

    label 'process_single'

    input:
    path overlap_tsv        // CATALOG_OVERLAP catalog_overlap.tsv
    path consensus_gtf      // TRANSCRIPT_LENGTH_FILTER kept GTF
    val bed_names           // catalog names (parallel to `beds`)
    path beds               // catalog BED files (may be empty)

    output:
    path "novelty_status.tsv", emit: novelty
    path "versions.yml", emit: versions

    script:
    def bed_args = bed_names.indices
        .collect { i -> "--bed ${bed_names[i]}=\"${beds[i]}\"" }
        .join(' ')
    """
    check_novelty.py \\
        --overlap "${overlap_tsv}" \\
        --consensus-gtf "${consensus_gtf}" \\
        ${bed_args} \\
        --out novelty_status.tsv

    emit_versions.sh "${task.process}" "${task.container}" python3 python3 --version
    """
}
