// F1.3a: per-sample featureCounts counts for the flank SAF.
// One task per BAM (both -s 2 and -s 1 passes), so the work parallelises across
// samples instead of serialising all BAMs in one featureCounts call.

process FLANK_COUNTS {
    tag "${sample}"

    input:
    tuple val(sample), path(bam)
    path flank_saf

    output:
    tuple val(sample), path("${sample}.flank_s2.tsv"), emit: s2
    tuple val(sample), path("${sample}.flank_s1.tsv"), emit: s1

    script:
    """
    featureCounts \\
        -F SAF \\
        -a "${flank_saf}" \\
        -o "${sample}.flank_s2.tsv" \\
        -T ${task.cpus} \\
        -p --countReadPairs \\
        -s 2 \\
        "${bam}"

    featureCounts \\
        -F SAF \\
        -a "${flank_saf}" \\
        -o "${sample}.flank_s1.tsv" \\
        -T ${task.cpus} \\
        -p --countReadPairs \\
        -s 1 \\
        "${bam}"
    """
}
