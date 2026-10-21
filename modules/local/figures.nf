// F5.1: catalog figures. Non-blocking (errorStrategy 'ignore'); each panel is
// also wrapped in tryCatch inside bin/figures.R so an odd table yields a
// manifest row instead of a failed run.

process FIGURES {

    label 'process_single'
    errorStrategy 'ignore'

    input:
    path catalog_genes      // deliverables catalog_genes.tsv
    path tier_by_class      // tier_final_by_class.tsv

    output:
    path "*.png", emit: figures
    path "figures_manifest.tsv", emit: manifest

    script:
    """
    figures.R --catalog "${catalog_genes}" --by-class "${tier_by_class}" --outdir .
    """
}
