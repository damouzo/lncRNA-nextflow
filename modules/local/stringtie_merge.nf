// A4: StringTie2 merge across all samples into single discovery GTF

process STRINGTIE2_MERGE {

    input:
    path(gtf_list)         // collected list of per-sample transcript GTFs
    path reference_gtf

    output:
    path "merged_transcripts.gtf", emit: merged_gtf
    path "mergelist.txt", emit: mergelist

    script:
    """
    printf '%s\\n' ${gtf_list.collect { "\"${it}\"" }.join(' ')} > mergelist.txt

    stringtie --merge \\
        -p ${task.cpus} \\
        -G "${reference_gtf}" \\
        -o merged_transcripts.gtf \\
        mergelist.txt
    """
}
