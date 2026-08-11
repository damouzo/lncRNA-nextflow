// Build CPAT species-specific models from the reference GTF
// Positive set: CDS from protein_coding genes
// Negative set: transcripts annotated as lncRNA
//
// Runs once per GTF version; output is cached for reuse by the CPAT module.
// Skipped if cpat_hexamer + cpat_logit_model are already provided in params.

process CPAT_BUILD_MODEL {
    label 'large_task'
    cpus = 8
    memory = 32.GB

    input:
    path genome_fa
    path reference_gtf

    output:
    path "cpat_models/Human_Hexamer.tsv",   emit: hexamer_table
    path "cpat_models/Human_logitModel.RData", emit: logit_model

    script:
    """
    mkdir -p cpat_models

    # --- Extract CDS from protein-coding genes (positive set) ---
    # Filter GTF to protein_coding genes, extract CDS
    gffread -g "${genome_fa}" \\
        -x cpat_models/coding_CDS.fa \\
        -C \\
        "${reference_gtf}"

    # --- Extract lncRNA transcript sequences (negative set) ---
    # Filter GTF to lncRNA biotype, extract full transcripts
    awk -F'\\t' '\$3 == "transcript"' "${reference_gtf}" \\
        | grep -E 'transcript_biotype "lncRNA"' \\
        > cpat_models/lncrna_transcripts.gtf

    if [ -s cpat_models/lncrna_transcripts.gtf ]; then
        gffread -w cpat_models/lncrna_transcripts.fa \\
            -g "${genome_fa}" \\
            cpat_models/lncrna_transcripts.gtf
    else
        echo "WARNING: no lncRNA transcripts found in GTF; using empty negative set"
        touch cpat_models/lncrna_transcripts.fa
    fi

    # --- Build CPAT hexamer frequency table ---
    make_hexamer_tab.py \\
        -c cpat_models/coding_CDS.fa \\
        -n cpat_models/lncrna_transcripts.fa \\
        > cpat_models/Human_Hexamer.tsv

    # --- Train logistic regression model ---
    make_logitModel.py \\
        -x cpat_models/Human_Hexamer.tsv \\
        -c cpat_models/coding_CDS.fa \\
        -n cpat_models/lncrna_transcripts.fa \\
        -o cpat_models/Human
    """
}
