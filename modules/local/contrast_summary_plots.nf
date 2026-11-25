// Cross-contrast comparison figures + tables, split by lncRNA subset
// (combined = annotated + novel / annotated / novel). Panels: significant counts
// per contrast (up/down), Venn (<= 4 contrasts) and UpSet of the significant
// lncRNAs, plus per-subset membership tables. Non-blocking: a missing plotting
// package or an empty input only skips panels, never fails the run.
process CONTRAST_SUMMARY_PLOTS {

    label 'process_single'
    errorStrategy 'ignore'

    input:
    path across          // lncrna_across_contrasts.tsv (per-contrast .padj blocks)
    path plot_script     // bin/contrast_summary_plots.R, staged so edits are hashed

    output:
    path "plots/**", optional: true
    path "tables/**", optional: true
    path "versions.yml", emit: versions

    script:
    """
    # `|| true` keeps a plotting hiccup from failing the task even with
    # errorStrategy 'ignore' (which would otherwise let a partial work dir through).
    Rscript contrast_summary_plots.R \\
        --across "${across}" \\
        --outdir . || true

    emit_versions.sh "${task.process}" "${task.container}" R R --version
    """
}
