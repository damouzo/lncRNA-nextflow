// DDX41/DHX34 positive-control report on the finished DE tables.
// Runs once, after every group's DESEQ2_DE, and writes run_info/qc_controls.tsv
// (folded into run_info by PUBLISH_RUN_INFO). Diagnostic only: it never fails
// the run, it just records where each control gene landed.

process QC_CONTROLS {

    input:
    path de_tables          // all per-contrast DE tables (staged flat)
    path de_manifest        // group<TAB>contrast<TAB>file, one row per DE table
    path contrasts_tsv      // contrast_name<TAB>numerator<TAB>denominator<TAB>batch
    path overview           // contrasts_overview.tsv (ordering: run after deliverables)

    output:
    path "qc_controls.tsv", emit: qc_controls
    path "versions.yml", emit: versions

    script:
    """
    qc_controls.py \\
        --manifest "${de_manifest}" \\
        --de-dir . \\
        --contrasts "${contrasts_tsv}" \\
        --outdir .

    emit_versions.sh "${task.process}" "${task.container}" python python3 --version
    """
}
