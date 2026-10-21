// F3.2: per-contrast design rank check (metadata only, before DESeq2).
// params.design_check = fail (default) | warn.

process DESIGN_CHECK {

    label 'process_single'

    input:
    tuple val(group), path(coldata), path(contrasts_tsv)
    val design_formula

    output:
    path "${group}_design_check.tsv", emit: design_check
    path "versions.yml", emit: versions

    script:
    """
    design_check.R \\
        --group "${group}" \\
        --coldata "${coldata}" \\
        --contrasts "${contrasts_tsv}" \\
        --design-formula "${design_formula}" \\
        --mode ${params.design_check} \\
        --outdir .

    mv design_check.tsv "${group}_design_check.tsv"
    emit_versions.sh "${task.process}" "${task.container}" R R --version
    """
}
