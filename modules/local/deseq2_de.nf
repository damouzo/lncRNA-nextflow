// DESeq2 differential expression.
// The R logic lives in bin/deseq2_contrasts.R (testable without Nextflow); this
// module is a thin wrapper. F3.5: the LFC cutoff is NOT passed to results() —
// it only affects the `significant` flag built in build_contrast_tables.py, on
// the shrunken log2FoldChange. F3.3: padj = NA diagnostics per contrast.
//
// Contrasts TSV: contrast_name<TAB>numerator<TAB>denominator<TAB>batch

process DESEQ2_DE {

    input:
    tuple val(group), path(counts_rds), path(coldata_csv), path(contrasts_tsv)
    val design_formula

    output:
    path "${group}_deseq2_results/*",                        emit: all_results
    path "${group}_deseq2_results/dds.rds",                  emit: dds_rds
    path "${group}_deseq2_results/design_diagnostics.txt",   emit: design_diagnostics
    tuple val(group), path("${group}_deseq2_results/padj_na_summary.tsv"), emit: padj_na_summary
    tuple val(group), path("${group}_deseq2_results/contrasts_summary.tsv"), emit: contrasts_summary
    tuple val(group), path("${group}_deseq2_results/*_DE.tsv"), emit: de_per_contrast
    path "versions.yml", emit: versions

    script:
    """
    deseq2_contrasts.R \\
        --group "${group}" \\
        --counts "${counts_rds}" \\
        --coldata "${coldata_csv}" \\
        --contrasts "${contrasts_tsv}" \\
        --design-formula "${design_formula}" \\
        --fdr ${params.fdr_threshold} \\
        --outdir .

    emit_versions.sh "${task.process}" "${task.container}" R R --version
    """
}
