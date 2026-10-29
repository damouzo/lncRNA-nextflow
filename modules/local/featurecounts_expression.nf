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

    // -O --fraction: candidate transcripts overlap heavily, so each read must be
    // split across every transcript it overlaps rather than stop at the first.
    // Deliberate, but it adds per-sample runtime and peak memory; revisit if the
    // candidate GTF grows.
    script:
    def strandness = Strandedness.featurecountsStrand(params.library_strandedness, params.featurecounts_strandness)
    """
    featureCounts \\
        -a "${candidate_gtf}" \\
        -o "${sample}.featureCounts.tsv" \\
        -t exon \\
        -g transcript_id \\
        -s ${strandness} \\
        -T ${task.cpus} \\
        -p \\
        --countReadPairs \\
        -O \\
        --fraction \\
        "${bam}"
    """
}
