// A2: Empirical strandedness detection via RSeQC infer_experiment.py
//
// Consumes a deterministic reference BED12 (BUILD_INFER_BED) so the `failed`
// fraction is not inflated by overlapping/antisense exons. The empirical result
// is checked against the declared contract by CHECK_STRANDEDNESS.

process INFER_STRANDEDNESS {
    tag "${sample}"

    input:
    tuple val(sample), path(bam)
    path ref_bed12

    output:
    tuple val(sample), path("${sample}.strandedness.txt"), emit: strandedness

    script:
    """
    infer_experiment.py \\
        -r "${ref_bed12}" \\
        -i "${bam}" \\
        -s 200000 \\
        -q 30 \\
        > "${sample}.strandedness.txt"
    """
}
