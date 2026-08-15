// A5: gffcompare classification against reference GTF
// Classes transcripts: u (intergenic), i (intronic), x (antisense), etc.

process GFFCOMPARE {

    input:
    path merged_gtf
    path reference_gtf

    output:
    path "gffcmp.*", emit: all_output
    path "gffcmp.annotated.gtf", emit: annotated_gtf
    path "gffcmp.loci", emit: loci
    path "gffcmp.tracking", emit: tracking
    path "gffcmp.class_codes.txt", emit: class_summary

    script:
    """
    gffcompare \\
        -r "${reference_gtf}" \\
        -o gffcmp \\
        "${merged_gtf}"

# Quick class-code counts for the log
    # class_code lives in the transcript attributes (col 9)
    grep 'class_code' gffcmp.annotated.gtf \
        | sed 's/.*class_code "//;s/".*//' \
        | sort | uniq -c | sort -rn \
        > gffcmp.class_codes.txt
    """
}
