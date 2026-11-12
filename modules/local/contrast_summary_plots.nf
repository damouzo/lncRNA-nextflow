// Cross-contrast summary figures: significant counts per contrast (up/down) and
// the UpSet overlap of lncRNAs significant across contrasts. Non-blocking: a
// missing plotting package or an empty input only skips panels, never fails the run.
process CONTRAST_SUMMARY_PLOTS {

    label 'process_single'
    errorStrategy 'ignore'

    input:
    path overview        // contrasts_overview.tsv (one row per contrast)
    path across          // lncrna_across_contrasts.tsv (per-contrast .padj blocks)
    path plot_script     // bin/contrast_summary_plots.R, staged so edits are hashed

    output:
    path "sig_counts_by_contrast.png", optional: true
    path "sig_counts_by_contrast.pdf", optional: true
    path "upset_lncrna.png", optional: true
    path "upset_lncrna.pdf", optional: true
    path "versions.yml", emit: versions

    script:
    """
    # `|| true` keeps a plotting hiccup from failing the task even with
    # errorStrategy 'ignore' (which would otherwise let a partial work dir through).
    Rscript contrast_summary_plots.R \\
        --overview "${overview}" \\
        --across "${across}" \\
        --outdir . || true

    emit_versions.sh "${task.process}" "${task.container}" R R --version
    """
}
