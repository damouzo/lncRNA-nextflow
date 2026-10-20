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
    def strand_flag = Strandedness.stringtieFlag(params.library_strandedness)
    // Build the command explicitly (no stripIndent ambiguity) so the reverse case
    // is byte-identical to the pre-contract command and Phase A -resume still hits.
    def flag_line = strand_flag ? "    ${strand_flag} \\\n" : ""
    "stringtie \\\n" +
    "    -p ${task.cpus} \\\n" +
    "    -G \"${reference_gtf}\" \\\n" +
    flag_line +
    "    -o \"${sample}.transcripts.gtf\" \\\n" +
    "    -A \"${sample}.gene_abund.tab\" \\\n" +
    "    -C \"${sample}.coverage.gtf\" \\\n" +
    "    \"${bam}\"\n"
}
