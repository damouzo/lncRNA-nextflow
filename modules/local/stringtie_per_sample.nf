// A3: Per-sample StringTie2 transcriptome assembly
// Genome-guided assembly using BAM files as input

process STRINGTIE2_PER_SAMPLE {
    tag "${sample}"

    input:
    tuple val(sample), path(bam)
    path reference_gtf

    output:
    tuple val(sample), path("${sample}.transcripts.gtf"), emit: transcript_gtf
    path "${sample}.gene_abund.tab", emit: gene_abundance
    path "${sample}.coverage.gtf", emit: coverage_gtf

    script:
    """
    stringtie \\
        -p ${task.cpus} \\
        -G "${reference_gtf}" \\
        -e \\
        -o "${sample}.transcripts.gtf" \\
        -A "${sample}.gene_abund.tab" \\
        -C "${sample}.coverage.gtf" \\
        "${bam}"
    """
}
