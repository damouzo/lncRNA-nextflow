// A2: Empirical strandedness detection via RSeQC infer_experiment.py
// Determines library strand orientation for downstream tools

process INFER_STRANDEDNESS {
    tag "${sample}"
    label 'small_task'
    cpus = 1
    memory = 4.GB

    input:
    tuple val(sample), path(bam)
    path reference_gtf

    output:
    tuple val(sample), path("${sample}.strandedness.txt"), emit: strandedness

    script:
    """
    infer_experiment.py \
        -r <(head -n 200000 "${reference_gtf}" | awk '\$3 == "exon"') \
        -i "${bam}" \
        -s ${task.cpus} \
        > "${sample}.strandedness.txt"
    """
}
