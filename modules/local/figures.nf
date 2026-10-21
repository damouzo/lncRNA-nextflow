// F5.1: catalog figures. Non-blocking (errorStrategy 'ignore'); each panel is
// also wrapped in tryCatch inside bin/figures.R so an odd table yields a
// manifest row instead of a failed run.
//
// Inputs beyond the catalog: the per-group TPM TSVs and the cross-group sample
// metadata (condition per sample) feed the F5.4 expression-by-position panel.
// Both are optional: with no TPM the panel is reported as skipped.

process FIGURES {

    label 'process_single'
    errorStrategy 'ignore'

    input:
    path catalog_genes      // deliverables catalog_genes.tsv (with lnc_position_class)
    path tier_by_class      // tier_final_by_class.tsv
    path tpm_tsvs           // collected deliverables/expression/<group>_tpm.tsv (may be empty)
    path sample_metadata    // deliverables/expression/sample_metadata.tsv (sample, condition)

    output:
    path "*.png", emit: figures
    path "figures_manifest.tsv", emit: manifest

    script:
    def tpm_args = tpm_tsvs.collect { "--tpm \"${it}\"" }.join(' ')
    """
    figures.R --catalog "${catalog_genes}" --by-class "${tier_by_class}" \
        ${tpm_args} --sample-metadata "${sample_metadata}" --outdir .
    """
}
