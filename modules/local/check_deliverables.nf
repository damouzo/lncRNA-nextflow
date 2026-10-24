// F7.3 guard: fail the run when the per-contrast stage did not receive every
// (group, contrast) that SPLIT_CONTRASTS declared. The ids are compared on the
// workflow channel values, not on the published tree: publishDir copies are
// asynchronous (dispatched to a thread pool after outputs are bound), so the
// delivered tree cannot be read reliably from within the same run. Depends on the
// overview so the master tables are also complete before the run ends.

process CHECK_DELIVERABLES {

    label 'process_single'

    input:
    val expected_contrasts    // declared <group>__<contrast> ids
    val delivered_contrasts   // ids that reached the per-contrast stage
    path overview             // contrasts_overview.tsv (ordering dependency)
    path contract             // docs/deliverables_contract.json

    output:
    path "versions.yml", emit: versions

    script:
    """
    validate_deliverables.py \\
        --dir "." \\
        --contract "${contract}" \\
        --expected-ids "${expected_contrasts.join(',')}" \\
        --found-contrasts "${delivered_contrasts.join(',')}" \\
        --contrast-coverage-only

    emit_versions.sh "${task.process}" "${task.container}" python3 python3 --version
    """
}
