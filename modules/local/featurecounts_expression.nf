// Quantify expression of candidate transcripts across all samples
// Uses featureCounts to count reads per transcript in the candidate GTF
// Produces a merged count matrix for expression recurrence filtering

process FEATURECOUNTS_EXPRESSION {

    input:
    tuple val(sample), path(bam)    // all validated BAMs (collected)
    path candidate_gtf

    output:
    path "${sample}.featureCounts.tsv", emit: per_sample_counts
    path "${sample}.featureCounts.tsv.summary", emit: summary

    script:
    """
    featureCounts \\
        -a "${candidate_gtf}" \\
        -o "${sample}.featureCounts.tsv" \\
        -t exon \\
        -g transcript_id \\
        -s 0 \\
        -T ${task.cpus} \\
        -p \\
        --countReadPairs \\
        "${bam}"
    """
}
