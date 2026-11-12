// F5.1: catalog figures. Non-blocking (errorStrategy 'ignore'); each panel is
// also wrapped in tryCatch inside bin/figures.R so an odd table yields a
// manifest row instead of a failed run.
//
// Inputs beyond the catalog: the per-group TPM TSVs and the cross-group sample
// metadata (condition per sample) feed the F5.4 expression-by-position panel;
// the coding consensus (CPAT/CPC2) and the candidate funnel feed the catalog
// characterization panels. All optional: a missing input is reported as skipped.

process FIGURES {

    label 'process_single'
    errorStrategy 'ignore'

    input:
    path catalog_genes      // deliverables catalog_genes.tsv (with lnc_position_class)
    path tier_by_class      // tier_final_by_class.tsv
    path tpm_tsvs           // collected deliverables/expression/<group>_tpm.tsv (may be empty)
    path sample_metadata    // run_info sample_index.tsv (sample, norm_group, condition)
    path comparisons        // the run's comparisons CSV (denominator = control condition)
    path coding_consensus   // CODING_CONSENSUS coding_consensus.tsv (CPAT/CPC2)
    path candidates_all     // deliverables/catalog candidates_all.tsv (discovery funnel)

    output:
    path "*.png", emit: figures
    path "*.pdf", emit: figures_pdf
    path "figures_manifest.tsv", emit: manifest

    script:
    def tpm_args = tpm_tsvs.collect { "--tpm \"${it}\"" }.join(' ')
    """
    figures.R --catalog "${catalog_genes}" --by-class "${tier_by_class}" \
        ${tpm_args} --sample-metadata "${sample_metadata}" \
        --comparisons "${comparisons}" \
        --coding-consensus "${coding_consensus}" \
        --candidates "${candidates_all}" \
        --cpat-threshold ${params.cpat_coding_threshold} \
        --cpc2-threshold ${params.cpc2_coding_threshold} \
        --outdir .
    """
}
