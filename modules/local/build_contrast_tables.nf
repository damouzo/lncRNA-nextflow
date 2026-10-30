// F4.2: per-contrast master tables. Last node of the graph.
// One row per gene of the quantified universe (significant or not, tested or not),
// so nothing is lost. The LFC cutoff is applied only to the `significant_padj_lfc1`
// flag here; it must never enter DESeq2 results() as lfcThreshold.

process BUILD_CONTRAST_TABLES {

    label 'process_single'

    input:
    tuple val(group), val(ct), path(de_tsv), path(universe_tsv)
    path catalog_genes      // deliverables/catalog catalog_genes.tsv (annotation)

    output:
    path "de_full.tsv"
    path "de_full_significant.tsv"
    path "de_lncrna_full.tsv"
    path "de_lncrna_significant.tsv"
    path "summary.tsv"
    // Keyed copy for the overview builder: the contrast id is not derivable from
    // the work-dir filenames, and the downstream task collects several of these.
    tuple val("${group}__${ct}"), path("summary.tsv"), path("de_full_significant.tsv"), emit: overview_input
    path "versions.yml", emit: versions

    script:
    """
    build_contrast_tables.py \\
        --universe "${universe_tsv}" \\
        --de "${de_tsv}" \\
        --catalog "${catalog_genes}" \\
        --contrast "${group}__${ct}" \\
        --fdr ${params.fdr_threshold} \\
        --lfc ${params.lfc_threshold} \\
        --extreme-basemean ${params.extreme_basemean} \\
        --genome-build "${params.genome_build ?: ''}" \\
        --out de_full.tsv \\
        --significant-out de_full_significant.tsv \\
        --lncrna-out de_lncrna_full.tsv \\
        --lncrna-significant-out de_lncrna_significant.tsv \\
        --summary-out summary.tsv

    emit_versions.sh "${task.process}" "${task.container}" python3 python3 --version
    """
}
