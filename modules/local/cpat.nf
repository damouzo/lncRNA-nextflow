// CPAT — Coding-Potential Assessment Tool
// Receives hexamer table + logit model (pre-built or auto-generated)

process CPAT {

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
        -d "${cpat_logit_model}" \\
        -x "${cpat_hexamer}" \\
        --antisense

    # Build ID/mRNA (yes/no) table consumed by CODING_CONSENSUS.
    # CPAT reports ORFs per transcript; a transcript is non-coding when it has
    # no ORF, or when its best-ORF coding probability is below the threshold.
    awk -F'\\t' -v thr=${params.cpat_coding_threshold} \
        'NR == 1 { next }
         { print \$1 "\\t" (\$11 < thr ? "no" : "yes") }' \
        cpat.ORF_prob.best.tsv > cpat_mrna.tsv
    awk '{ print \$1 "\\tno" }' cpat.no_ORF.txt >> cpat_mrna.tsv
    sort -u cpat_mrna.tsv > cpat_results.tmp
    printf 'ID\\tmRNA\\n' > cpat_results.tsv
    cat cpat_results.tmp >> cpat_results.tsv
    """
}
