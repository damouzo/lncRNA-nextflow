// A8: CPC2 — Coding Potential Calculator 2
// Score < 0.5 = non-coding

process CPC2 {

    input:
    path length_filtered_gtf
    path genome_fa

    output:
    path "cpc2_results.txt", emit: cpc2_output

    script:
    """
    gffread -w cpc2_transcripts.fa -g "${genome_fa}" "${length_filtered_gtf}"

    python3 `command -v CPC2.py` \\
        -i cpc2_transcripts.fa \\
        -o cpc2_results
    """
}
