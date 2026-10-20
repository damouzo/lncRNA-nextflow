// F1.5b + F1.6: transcript-level tier_final, gene aggregation and sensitivity.
// Reads the metrics + derived thresholds; terminal (nothing downstream).

process COMPUTE_TIER_FINAL {

    input:
    path frozen_gtf
    path reference_gtf
    path flank_metrics
    path antisense_metrics
    path thresholds
    path sensitivity

    output:
    path "tier_final_transcripts.tsv", emit: tier_final_transcripts
    path "tier_final_genes.tsv", emit: tier_final_genes
    path "tier_final_sensitivity.tsv", emit: sensitivity_table

    script:
    def allow_single_exon = params.tier1_allow_single_exon ? "--tier1-allow-single-exon" : ""
    """
    f1_tier_final.py \\
        --frozen-gtf "${frozen_gtf}" \\
        --reference-gtf "${reference_gtf}" \\
        --flank-metrics "${flank_metrics}" \\
        --antisense-metrics "${antisense_metrics}" \\
        --thresholds "${thresholds}" \\
        --sensitivity-thresholds "${sensitivity}" \\
        --outdir . \\
        --tier1-classes '${params.tier1_classes.join(",")}' \\
        --tier1a-classes '${params.tier1a_classes.join(",")}' \\
        --tier1-parent-biotypes '${params.tier1_parent_biotypes.join(",")}' \\
        ${allow_single_exon}
    """
}
