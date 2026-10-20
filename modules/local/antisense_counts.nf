// F1.4a: per-sample featureCounts counts for the antisense overlap SAF.
// One task per BAM (both -s 2 and -s 1 passes).

process ANTISENSE_COUNTS {
    tag "${sample}"

    input:
    tuple val(sample), path(bam)
    path antisense_saf

    output:
    tuple val(sample), path("${sample}.anti_s2.tsv"), emit: s2
    tuple val(sample), path("${sample}.anti_s1.tsv"), emit: s1

    script:
    """
    featureCounts \\
        -F SAF \\
        -a "${antisense_saf}" \\
        -o "${sample}.anti_s2.tsv" \\
        -T ${task.cpus} \\
        -p --countReadPairs \\
        -s 2 \\
        "${bam}"

    featureCounts \\
        -F SAF \\
        -a "${antisense_saf}" \\
        -o "${sample}.anti_s1.tsv" \\
        -T ${task.cpus} \\
        -p --countReadPairs \\
        -s 1 \\
        "${bam}"
    """
}
