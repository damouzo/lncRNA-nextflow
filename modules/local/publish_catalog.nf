// F4.3 / F4.3b: publish the catalog deliverable.
// Builds the contract tables under deliverables/catalog/ (catalog_genes, catalog_tx,
// candidates_all) with bin/build_catalog_tables.py, and the annotation/ (GTF with LNC
// ids + FASTA via gffread). No BED12 (intentionally not published).
// Tier tables (tier_final_by_class/sensitivity) are published by COMPUTE_TIER_FINAL.

process PUBLISH_CATALOG {

    label 'process_single'

    input:
    path gene_catalog           // BUILD_GENE_CATALOG (full universe)
    path tier_final_genes       // COMPUTE_TIER_FINAL, gene level
    path tier_final_transcripts // COMPUTE_TIER_FINAL, transcript level
    path gffcmp_gtf             // GFFCOMPARE annotated GTF (candidate funnel)
    path recurrence             // EXPRESSION_RECURRENCE recurrence_filtered.tsv
    path frozen_gtf             // ANNOTATION_FREEZE frozen_lncrna.gtf (frozen id set)
    path catalog_gtf            // ASSIGN_CATALOG_IDS frozen_lncrna.catalog.gtf (LNC ids)
    path genome_fa

    output:
    path "catalog_genes.tsv"
    path "catalog_tx.tsv"
    path "candidates_all.tsv"
    path "annotation/lncrna_catalog.gtf"
    path "annotation/lncrna_catalog.fa"

    script:
    """
    build_catalog_tables.py \\
        --gene-catalog "${gene_catalog}" \\
        --tier-final-genes "${tier_final_genes}" \\
        --tier-final-transcripts "${tier_final_transcripts}" \\
        --gffcmp-gtf "${gffcmp_gtf}" \\
        --recurrence "${recurrence}" \\
        --frozen-gtf "${frozen_gtf}" \\
        --min-length ${params.min_transcript_length} \\
        --outdir .

    mkdir -p annotation
    cp "${catalog_gtf}" annotation/lncrna_catalog.gtf
    gffread -w annotation/lncrna_catalog.fa -g "${genome_fa}" annotation/lncrna_catalog.gtf
    """
}
