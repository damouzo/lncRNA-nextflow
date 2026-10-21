// B5 / F3.1: cis-regulatory candidate associations (novel + annotated lncRNAs).
// The R logic lives in bin/cis_associations.R (testable without Nextflow); this
// module is a thin wrapper. params.cis_normalize selects the normalization
// (VST + removeBatchEffect vs log2(counts + 1)).

process CIS_ASSOCIATIONS {

    input:
    tuple val(group), path(counts_rds), path(deseq2_summary), path(coldata_csv)
    path gene_catalog       // gene_catalog.tsv: origin + coordinates for every gene

    output:
    tuple val(group), path("${group}_cis_pairs.tsv"),              emit: cis_pairs
    tuple val(group), path("${group}_cis_pairs_significant.tsv"),  emit: cis_pairs_sig
    path "${group}_cis_associations_summary.tsv", emit: summary
    path "versions.yml", emit: versions

    script:
    """
    cis_associations.R \\
        --group "${group}" \\
        --gene-catalog "${gene_catalog}" \\
        --counts "${counts_rds}" \\
        --coldata "${coldata_csv}" \\
        --cis-window ${params.cis_window} \\
        --fdr ${params.fdr_threshold} \\
        --normalize ${params.cis_normalize} \\
        --outdir .

    emit_versions.sh "${task.process}" "${task.container}" R R --version
    """
}
