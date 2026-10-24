// F4.2/F4.3: single-task overview of every contrast (concatenation of the
// per-contrast summaries plus the novel Tier-1 significant count). The inputs are
// already renamed by contrast id in the workflow, so they can be staged together.

process BUILD_CONTRASTS_OVERVIEW {

    label 'process_single'

    input:
    path summaries        // collected <id>__summary.tsv
    path significants     // collected <id>__significant.tsv
    path catalog_genes

    output:
    path "contrasts_overview.tsv", emit: overview
    path "versions.yml", emit: versions

    script:
    """
    build_contrasts_overview.py \\
        --summaries-dir . \\
        --catalog "${catalog_genes}" \\
        --out contrasts_overview.tsv

    emit_versions.sh "${task.process}" "${task.container}" python3 python3 --version
    """
}
