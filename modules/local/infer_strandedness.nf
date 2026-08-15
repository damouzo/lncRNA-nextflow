// A2: Empirical strandedness detection via RSeQC infer_experiment.py
// Determines library strand orientation for downstream tools

process INFER_STRANDEDNESS {
    tag "${sample}"

    input:
    tuple val(sample), path(bam)
    path reference_gtf

    output:
    tuple val(sample), path("${sample}.strandedness.txt"), emit: strandedness

script:
    """
    gffread --bed "${reference_gtf}" -o ref.bed 2>/dev/null
    awk -F'\\t' 'NF >= 12 {print \$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$5"\\t"\$6"\\t"\$7"\\t"\$8"\\t"\$9"\\t"\$10"\\t"\$11"\\t"\$12}' ref.bed > ref_exon.bed12

    infer_experiment.py \\
        -r ref_exon.bed12 \\
        -i "${bam}" \\
        -s 200000 \\
        -q 30 \\
        > "${sample}.strandedness.txt"
    """
}
