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
    head -n 200000 "${reference_gtf}" | awk '\$3 == "exon"' > exon_ref.bed
    infer_experiment.py \\
        -r exon_ref.bed \\
        -i "${bam}" \\
        -s ${task.cpus} \\
        > "${sample}.strandedness.txt"
    """
}
