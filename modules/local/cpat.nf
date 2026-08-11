// CPAT — Coding-Potential Assessment Tool
// Receives hexamer table + logit model (pre-built or auto-generated)

process CPAT {
    label 'medium_task'

    input:
    path length_filtered_gtf
    path genome_fa
    tuple path(cpat_hexamer), path(cpat_logit_model)

    output:
    path "cpat_results.tsv", emit: cpat_output

    script:
    """
    gffread -w cpat_transcripts.fa -g "${genome_fa}" "${length_filtered_gtf}"

    cpat.py \\
        -g cpat_transcripts.fa \\
        -o cpat \\
        -d "${cpat_hexamer}" \\
        -x "${cpat_logit_model}" \\
        --antisense

    awk -F'\\t' 'NR==1 || \$6 == "no"' cpat.ORF_prob.tsv > cpat_results.tsv
    """
}
