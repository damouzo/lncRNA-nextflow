// Expression recurrence filter (condition-blind)
// Reads featureCounts matrices from all samples, computes per-transcript
// expression prevalence, and keeps transcripts expressed in >= N replicates.
// The CPM denominator is the total fragments per sample (featureCounts .summary,
// all categories), not the candidate-only column sums (H2 / task 3.1).

process EXPRESSION_RECURRENCE {

    input:
    path catalog_table        // from catalog_overlap
    path fc_tsv               // collected featureCounts TSVs from all samples
    path fc_summaries         // collected featureCounts .summary files (total fragments)
    path merged_gtf           // merged StringTie2 GTF (for pass-through)

    output:
    path "recurrence_filtered.tsv", emit: recurrence_table
    path "recurrence_filtered.gtf", emit: recurrence_gtf
    path "recurrence_stats.txt", emit: stats

    script:
    """
    expression_recurrence.R \\
        --catalog "${catalog_table}" \\
        --merged-gtf "${merged_gtf}" \\
        --fc-dir . \\
        --min-reps ${params.expression_recurrence_min_reps} \\
        --cpm-threshold ${params.recurrence_cpm_threshold} \\
        --outdir .
    """
}
