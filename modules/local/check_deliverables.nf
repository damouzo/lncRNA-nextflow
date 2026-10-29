// F7.3 guard: fail the run when the per-contrast stage did not receive every
// (group, contrast) that SPLIT_CONTRASTS declared. The ids are compared on the
// workflow channel values, not on the published tree: publishDir copies are
// asynchronous (dispatched to a thread pool after outputs are bound), so the
// delivered tree cannot be read reliably from within the same run. Depends on the
// overview so the master tables are also complete before the run ends.

process CHECK_DELIVERABLES {

    label 'process_single'

    input:
    val expected_contrasts    // declared <group>__<contrast> ids (applicability-resolved)
    val delivered_contrasts   // ids that reached the per-contrast stage
    path declared_contrasts   // the raw contrasts file: count of declared rows
    path overview             // contrasts_overview.tsv (ordering dependency)
    path contract             // docs/deliverables_contract.json

    output:
    path "versions.yml", emit: versions

    script:
    """
    # 0.2 count guard: the number of delivered contrasts must equal the number of
    # rows declared in the contrasts file. This is computed from the raw file, not
    # from the SPLIT_CONTRASTS output, so a contrast silently dropped upstream is
    # still required to be delivered.
    validate_deliverables.py \\
        --dir "." \\
        --contract "${contract}" \\
        --expect-contrasts "${declared_contrasts}" \\
        --found-contrasts "${delivered_contrasts.join(',')}" \\
        --contrast-coverage-only

    # id guard: every applicability-resolved id must have been delivered.
    validate_deliverables.py \\
        --dir "." \\
        --contract "${contract}" \\
        --expected-ids "${expected_contrasts.join(',')}" \\
        --found-contrasts "${delivered_contrasts.join(',')}" \\
        --contrast-coverage-only

    emit_versions.sh "${task.process}" "${task.container}" python3 python3 --version
    """
}
