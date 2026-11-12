// F4.3 / F4.3b: publish the catalog deliverable.
// Builds the contract tables under deliverables/catalog/ (catalog_genes, catalog_tx,
// candidates_all) with bin/build_catalog_tables.py, and the annotation/ (GTF with LNC
// ids + FASTA via gffread + a short README). No BED12 (intentionally not published).
// Tier tables (tier_final_by_class/sensitivity) are published by COMPUTE_TIER_FINAL.

process PUBLISH_CATALOG {

    label 'process_single'

    input:
    path gene_catalog           // BUILD_GENE_CATALOG (full universe)
    path tier_final_genes       // COMPUTE_TIER_FINAL, gene level
    path tier_final_transcripts // COMPUTE_TIER_FINAL, transcript level
    path gffcmp_gtf             // GFFCOMPARE annotated GTF (candidate funnel)
    path length_filtered_gtf    // TRANSCRIPT_LENGTH_FILTER kept GTF (parity check)
    path recurrence             // EXPRESSION_RECURRENCE recurrence_filtered.tsv
    path frozen_gtf             // ANNOTATION_FREEZE frozen_lncrna.gtf (frozen id set)
    path catalog_gtf            // ASSIGN_CATALOG_IDS frozen_lncrna.catalog.gtf (LNC ids)
    path genome_fa

    output:
    path "catalog_genes.tsv", emit: catalog_genes
    path "catalog_tx.tsv"
    path "candidates_all.tsv"
    path "novel_lncrnas.tsv"
    path "annotation/lncrna_catalog.gtf"
    path "annotation/lncrna_catalog.fa"
    path "annotation/README.md"
    path "versions.yml", emit: versions

    script:
    """
    build_catalog_tables.py \\
        --gene-catalog "${gene_catalog}" \\
        --tier-final-genes "${tier_final_genes}" \\
        --tier-final-transcripts "${tier_final_transcripts}" \\
        --catalog-gtf "${catalog_gtf}" \\
        --gffcmp-gtf "${gffcmp_gtf}" \\
        --length-filtered-gtf "${length_filtered_gtf}" \\
        --recurrence "${recurrence}" \\
        --frozen-gtf "${frozen_gtf}" \\
        --min-length ${params.min_transcript_length} \\
        --novel-classes '${params.novel_classes.join(",")}' \\
        --isoform-classes '${params.isoform_classes.join(",")}' \\
        --outdir .

    mkdir -p annotation
    cp "${catalog_gtf}" annotation/lncrna_catalog.gtf
    gffread -w annotation/lncrna_catalog.fa -g "${genome_fa}" annotation/lncrna_catalog.gtf

    cat > annotation/README.md <<'EOF'
# annotation/

Frozen lncRNA catalog in annotation format.

- `lncrna_catalog.gtf` - the frozen novel lncRNAs with stable ids. Gene and
  transcript ids are `LNC000001` / `LNC000001.1`; each record also carries
  `catalog_id`, `mstrg_gene_id` / `mstrg_transcript_id` (the discovery ids)
  and `class_code` / `cmp_ref` from gffcompare.
- `lncrna_catalog.fa` - transcript sequences for the same set (gffread -w);
  headers match the `catalog_transcript_id` of `../catalog_tx.tsv`.

Novel features only (no reference genes). Coordinate build follows the run's
reference. Not a filter result: the identity of every catalog transcript.
EOF

    emit_versions.sh "${task.process}" "${task.container}" gffread gffread --version
    """
}
