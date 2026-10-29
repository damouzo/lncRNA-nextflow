// B3: tximport + independent filtering (condition-blind).
// The R logic now lives in bin/tximport_filter.R (testable without Nextflow); this
// module is a thin wrapper. Same core outputs as before plus the TSV expression
// tables and the pre-filter universe with a reason per gene (F4.4 / 5f).

process TXIMPORT_AND_FILTER {

    input:
    tuple val(group), path(quant_dirs), path(coldata_csv)  // keyed: no positional pairing downstream
    path tx2gene_detailed
    path catalog_genes    // catalog_genes.tsv: gene_group + catalog_id for the TSV tables

    output:
    tuple val(group), path("${group}_tximport_counts.rds"), emit: counts_rds
    // 2.1: tximport-derived DESeqDataSet (kept genes) carrying the length offset.
    tuple val(group), path("${group}_tximport_dds.rds"),    emit: dds_rds
    // Keyed by group: downstream consumers .join() it per norm_group instead of
    // pairing single-execution/queue channels positionally across >1 norm_group.
    tuple val(group), path("${group}_tximport_tpm.rds"),    emit: tpm_rds
    // H0: coldata aligned to the counts columns (samplesheet order != quant-dir
    // order). Consumers must use this instead of SPLIT_BY_NORMGROUP's coldata.
    // Named coldata_aligned_* so it never collides with the symlinked input
    // coldata_<group>.csv that Nextflow stages into this task dir.
    tuple val(group), path("coldata_aligned_${group}.csv"), emit: aligned_coldata
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
        --catalog "${catalog_genes}" \\
        --min-count ${params.tximport_min_count} \\
        --outdir .

    emit_versions.sh "${task.process}" "${task.container}" R R --version
    emit_r_versions.sh "${task.process}" "${task.container}" DESeq2 tximport
    """
}
