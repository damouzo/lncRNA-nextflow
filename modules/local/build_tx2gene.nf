// BUILD_TX2GENE — the single producer of tx2gene_detailed.tsv for the whole run.
//
// The quantification universe (quant_universe.gtf) holds the catalog transcripts
// AND the reference-derived novel isoforms (isoform_classes, e.g. `j`). This
// module builds the canonical transcript->gene map the DE universe is defined
// from:
//   * reference transcripts      -> their reference gene;
//   * novel catalog transcripts  -> their own MSTRG gene;
//   * novel isoform transcripts  -> the REFERENCE PARENT gene (ref_gene_id).
//
// Remapping the isoforms to the parent is what keeps the parent gene's counts
// whole: transcript-level EM still splits reads between the reference isoforms
// and the novel isoform, and tximport sums them back under one gene_id. Without
// it the parent is diluted by a second, near-identical MSTRG gene.
//
// No other process may construct tx2gene_detailed.tsv (enforced by
// bin/tests/check_wiring.py). Isoforms are reported in novel_isoforms.tsv, never
// as novel lncRNAs, and never registered in the catalog.

process BUILD_TX2GENE {

    label 'process_single'

    input:
    path quant_universe_gtf    // ANNOTATION_FREEZE quant_universe.gtf (catalog + isoforms)
    path reference_gtf

    output:
    path "tx2gene_detailed.tsv", emit: tx2gene_detailed
    path "novel_isoforms.tsv",   emit: novel_isoforms
    path "tx2gene_manifest.json", emit: manifest
    path "versions.yml", emit: versions

    script:
    """
    build_tx2gene.R \\
        --reference "${reference_gtf}" \\
        --universe "${quant_universe_gtf}" \\
        --isoform-classes "${params.isoform_classes.join(',')}" \\
        --outdir .

    emit_versions.sh "${task.process}" "${task.container}" R R --version
    emit_r_versions.sh "${task.process}" "${task.container}" readr dplyr
    """
}
