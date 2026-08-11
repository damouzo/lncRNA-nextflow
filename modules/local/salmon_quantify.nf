// B2: Salmon quantification with decoy-aware index
// Uses empirically determined strand parameter from Phase A

process SALMON_QUANTIFY {
    tag "${sample}"
    label 'large_task'
    cpus = 8
    memory = 32.GB

    input:
    tuple val(sample), val(condition), val(batch), path(fastq_1), path(fastq_2), path(bam), val(extra)
    path salmon_index

    output:
    tuple val(sample),
          path("${sample}_quant"),
          emit: quant_results

    script:
    def strand_flag = params.salmon_libtype ?: 'A'
    """
    salmon quant \\
        -i salmon_index \\
        -l ${strand_flag} \\
        -1 "${fastq_1}" \\
        -2 "${fastq_2}" \\
        --gcBias \\
        --seqBias \\
        --validateMappings \\
        -p ${task.cpus} \\
        -o "${sample}_quant"

    # Rename quant files for clarity
    cp "${sample}_quant/quant.sf" "${sample}_quant/quant.sf.bak" 2>/dev/null || true
    """
}
