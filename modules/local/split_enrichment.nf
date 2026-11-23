// Repartition the group-level enrichment table into per-contrast deliverables.
// One task per (group, contrast): the input keys are already that pair, so the
// publishDir can land the files in contrasts/<group>__<ct>/ from the task inputs.

process SPLIT_ENRICHMENT {

    label 'process_single'

    input:
    tuple val(group), val(ct), path(ora)

    output:
    path "enrichment_ora.tsv"
    path "versions.yml", emit: versions

    script:
    """
    split_enrichment.py \\
        --ora "${ora}" \\
        --contrast "${ct}" \\
        --outdir .

    emit_versions.sh "${task.process}" "${task.container}" python3 python3 --version
    """
}
