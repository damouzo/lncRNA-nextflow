// DESeq2 differential expression.
// The R logic lives in bin/deseq2_contrasts.R (testable without Nextflow); this
// module is a thin wrapper. F3.5: the LFC cutoff is NOT passed to results() —
// it only affects the `significant` flag built in build_contrast_tables.py, on
// the shrunken log2FoldChange. F3.3: padj = NA diagnostics per contrast.
//
// Contrasts TSV: contrast_name<TAB>numerator<TAB>denominator<TAB>batch

process DESEQ2_DE {

    input:
    tuple val(group), path(counts_rds), path(tximport_dds), path(coldata_csv), path(contrasts_tsv)
    val design_formula

    output:
    path "${group}_deseq2_results/*",                        emit: all_results
    path "${group}_deseq2_results/dds.rds",                  emit: dds_rds
    path "${group}_deseq2_results/design_diagnostics.txt",   emit: design_diagnostics
    tuple val(group), path("${group}_deseq2_results/padj_na_summary.tsv"), emit: padj_na_summary
    tuple val(group), path("${group}_deseq2_results/contrasts_summary.tsv"), emit: contrasts_summary
    tuple val(group), path("${group}_deseq2_results/*_DE.tsv"), emit: de_per_contrast
    // 2.1: per-group flag saying whether the tximport length offset was used.
    // Group-prefixed by the script so the collected basenames stay unique.
    tuple val(group), path("${group}_deseq2_results/${group}_length_offset.txt"), emit: length_offset
    path "versions.yml", emit: versions

    script:
    """
    deseq2_contrasts.R \\
        --group "${group}" \\
        --counts "${counts_rds}" \\
        --dds "${tximport_dds}" \\
        --coldata "${coldata_csv}" \\
        --contrasts "${contrasts_tsv}" \\
        --design-formula "${design_formula}" \\
        --fdr ${params.fdr_threshold} \\
        --allow-no-offset ${params.allow_no_offset} \\
        --outdir .

    emit_versions.sh "${task.process}" "${task.container}" R R --version
    emit_r_versions.sh "${task.process}" "${task.container}" DESeq2 apeglm ashr
    """
}
