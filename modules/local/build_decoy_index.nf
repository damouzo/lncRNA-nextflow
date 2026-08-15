// Build decoy-aware Salmon index combining frozen lncRNAs + reference transcriptome
// If gencode_transcripts_fa is not provided, derives the reference transcriptome
// from the GTF annotation + genome FASTA via gffread.

process BUILD_DECOY_INDEX {

    input:
    path frozen_gtf
    path genome_fa
    val  gencode_transcripts_fa

    output:
    path "salmon_index/*", emit: index_dir

    script:
    """
    mkdir -p salmon_index

    # Extract novel lncRNA sequences from frozen GTF
    gffread -w lncrna_transcripts.fa -g "${genome_fa}" "${frozen_gtf}"

    # Determine reference transcriptome source
    if [ -n "\${gencode_transcripts_fa}" ] && [ -f "\${gencode_transcripts_fa}" ]; then
        ln -s "\${gencode_transcripts_fa}" ref_transcripts.fa
    else
        echo "Deriving reference transcriptome from GTF + genome..."
        gffread -w ref_transcripts.fa -g "${genome_fa}" "${params.gtf}"
    fi

    cat ref_transcripts.fa lncrna_transcripts.fa > combined_transcripts.fa

    grep "^>" "${genome_fa}" | cut -d ' ' -f 1 | sed 's/>//' > decoys.txt
    cat combined_transcripts.fa "${genome_fa}" > gentrome.fa

    salmon index \\
        -t gentrome.fa \\
        -d decoys.txt \\
        -p ${task.cpus} \\
        -i salmon_index \\
        --gencode
    """
}
