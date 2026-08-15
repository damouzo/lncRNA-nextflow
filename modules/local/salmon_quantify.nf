// B2: Salmon quantification with decoy-aware index
// Supports multi-lane FASTQs: fastq_1 and fastq_2 can be lists of Path objects
// Salmon's -1/-2 natively accept space-separated file lists — no concatenation needed

process SALMON_QUANTIFY {
    tag "${sample}"

    input:
    // fastq_1 and fastq_2 can be single files or lists (multi-lane)
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
        -1 ${fastq_1} \\
        -2 ${fastq_2} \\
        --gcBias \\
        --seqBias \\
        --validateMappings \\
        -p ${task.cpus} \\
        -o "${sample}_quant"

    # Rename quant files for clarity
    cp "${sample}_quant/quant.sf" "${sample}_quant/quant.sf.bak" 2>/dev/null || true
    """
}