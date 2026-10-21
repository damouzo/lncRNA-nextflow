// Slice D: one single-contrast file per (norm_group, contrast) that applies to
// the group. Lets REPORT_CONTRAST (and downstream builders) work per contrast
// with a keyed channel, and lets publishDir target contrasts/<grupo>__<contraste>/
// from the task inputs. Also rejects contrast names containing "__" or "/".

process SPLIT_CONTRASTS {

    label 'process_single'

    input:
    tuple val(group), path(coldata), path(contrasts_tsv)

    output:
    tuple val(group), path("contrasts/${group}__*.tsv"), emit: per_contrast

    script:
    """
    split_contrasts.py \\
        --coldata "${coldata}" \\
        --contrasts "${contrasts_tsv}" \\
        --group "${group}" \\
        --outdir contrasts
    """
}
