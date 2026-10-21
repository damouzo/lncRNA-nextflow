// B3: tximport + independent filtering (condition-blind).
// The R logic now lives in bin/tximport_filter.R (testable without Nextflow); this
// module is a thin wrapper. Same core outputs as before plus the TSV expression
// tables and the pre-filter universe with a reason per gene (F4.4 / 5f).

process TXIMPORT_AND_FILTER {

    input:
    tuple val(group), path(quant_dirs), path(coldata_csv)  // keyed: no positional pairing downstream
    path tx2gene_detailed

    output:
    tuple val(group), path("${group}_tximport_counts.rds"), emit: counts_rds
    // Keyed by group: downstream consumers .join() it per norm_group instead of
    // pairing single-execution/queue channels positionally across >1 norm_group.
    tuple val(group), path("${group}_tximport_tpm.rds"),    emit: tpm_rds
    path "${group}_sample_qc_pca.png",                      emit: pca_plot
    path "${group}_filter_stats.tsv",                       emit: filter_stats
    tuple val(group), path("${group}_counts.tsv"),          emit: counts_tsv
    tuple val(group), path("${group}_tpm.tsv"),             emit: tpm_tsv
    tuple val(group), path("${group}_universe.tsv"),        emit: universe_tsv
    path "versions.yml",                                    emit: versions

    script:
    """
    tximport_filter.R \\
        --group "${group}" \\
        --coldata "${coldata_csv}" \\
        --tx2gene "${tx2gene_detailed}" \\
        --outdir .

    emit_versions.sh "${task.process}" "${task.container}" R R --version
    """
}
