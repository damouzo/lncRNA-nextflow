// Repartition the group-level enrichment tables into per-contrast deliverables.
// One task per (group, contrast): the input keys are already that pair, so the
// publishDir can land the files in contrasts/<group>__<ct>/ from the task inputs.

process SPLIT_ENRICHMENT {

    label 'process_single'

    input:
    tuple val(group), val(ct), path(ora), path(gsea)

    output:
    path "enrichment_ora.tsv"
    path "enrichment_gsea.tsv"
    path "versions.yml", emit: versions

    script:
    """
    split_enrichment.py \\
        --ora "${ora}" \\
        --gsea "${gsea}" \\
        --contrast "${ct}" \\
        --outdir .

    emit_versions.sh "${task.process}" "${task.container}" python3 python3 --version
    """
}
