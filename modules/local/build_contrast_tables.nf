// F4.2: per-contrast master tables. Last node of the graph.
// One row per gene of the quantified universe (significant or not, tested or not),
// so nothing is lost. The LFC cutoff is applied only to the `significant_padj_lfc1`
// flag here; it must never enter DESeq2 results() as lfcThreshold.

process BUILD_CONTRAST_TABLES {

    label 'process_single'

    input:
    tuple val(group), val(ct), path(de_tsv), path(universe_tsv), path(tpm_tsv), path(coldata_csv)
    path catalog_genes      // deliverables/catalog catalog_genes.tsv (annotation)
    // bin/build_contrast_tables.py staged as an input so its content is part of the
    // task hash: Nextflow does not track bin/, so without this an edit to the script
    // would be silently served from cache under -resume.
    path build_script

    output:
    path "de_full.tsv"
    path "de_full_significant.tsv"
    path "de_lncrna_full.tsv"
    path "de_lncrna_significant.tsv"
    path "de_lncrna_significant_lfc1.tsv"
    path "summary.tsv"
    // Keyed copy for the overview builder: the contrast id is not derivable from
    // the work-dir filenames, and the downstream task collects several of these.
    tuple val("${group}__${ct}"), path("summary.tsv"), path("de_full_significant.tsv"), emit: overview_input
    tuple val("${group}__${ct}"), path("de_lncrna_full.tsv"), emit: lncrna_full
    path "versions.yml", emit: versions

    script:
    """
    # The contrast id is <numerator>_vs_<denominator> (e.g. MSC_DDX41_vs_Healthy).
    # Build it as a shell variable and split it in bash: no per-contrast file is
    # needed and no contrast can be dropped by a join further upstream.
    CID="${group}__${ct}"
    NUM="\${CID%%_vs_*}"
    DEN="\${CID#*_vs_}"

    build_contrast_tables.py \\
        --universe "${universe_tsv}" \\
        --de "${de_tsv}" \\
        --catalog "${catalog_genes}" \\
        --contrast "\$CID" \\
        --fdr ${params.fdr_threshold} \\
        --lfc ${params.lfc_threshold} \\
        --extreme-basemean ${params.extreme_basemean} \\
        --genome-build "${params.genome_build ?: ''}" \\
        --tpm "${tpm_tsv}" \\
        --coldata "${coldata_csv}" \\
        --numerator "\$NUM" \\
        --denominator "\$DEN" \\
        --out de_full.tsv \\
        --significant-out de_full_significant.tsv \\
        --lncrna-out de_lncrna_full.tsv \\
        --lncrna-significant-out de_lncrna_significant.tsv \\
        --lncrna-significant-lfc1-out de_lncrna_significant_lfc1.tsv \\
        --summary-out summary.tsv

    emit_versions.sh "${task.process}" "${task.container}" python3 python3 --version
    """
}
