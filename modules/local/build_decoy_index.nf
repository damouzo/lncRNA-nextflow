// Build decoy-aware Salmon index combining frozen lncRNAs + reference transcriptome
// If gencode_transcripts_fa is not provided, derives the reference transcriptome
// from the GTF annotation + genome FASTA via gffread.

process BUILD_DECOY_INDEX {

    input:
    path frozen_gtf
    path genome_fa
    val  gencode_transcripts_fa

    output:
    path "salmon_index", emit: index_dir

    script:
    """
    mkdir -p salmon_index

    # Extract novel lncRNA sequences from frozen GTF
    gffread -w lncrna_transcripts.fa -g "${genome_fa}" "${frozen_gtf}"

    # Determine reference transcriptome source
    # gencode_transcripts_fa is a `val` input, so it must be interpolated here
    REF_FA="${gencode_transcripts_fa}"
    if [ -n "\$REF_FA" ] && [ -f "\$REF_FA" ]; then
        ln -s "\$REF_FA" ref_transcripts.fa
    else
        echo "Deriving reference transcriptome from GTF + genome..."
        gffread -w ref_transcripts.fa -g "${genome_fa}" "${params.gtf}"
    fi

# Drop novel lncRNAs already present in the reference transcriptome
    # (Salmon rejects duplicate headers). IDs compared without version suffix.
    awk 'NR==FNR { if (\$1 ~ /^>/) { split(substr(\$1,2), a, "."); ref[a[1]] = 1 }; next }
         /^>/{ split(substr(\$1,2), a, "."); keep = !(a[1] in ref) }
         keep' ref_transcripts.fa lncrna_transcripts.fa > novel_only.fa

    grep "^>" "${genome_fa}" | cut -d ' ' -f 1 | sed 's/>//' > decoys.txt
    # Decoy-aware gentrome: full reference + novel lncRNAs + genome decoys
    cat ref_transcripts.fa novel_only.fa "${genome_fa}" > gentrome.fa

    salmon index \\
        -t gentrome.fa \\
        -d decoys.txt \\
        -p ${task.cpus} \\
        ${params.salmon_index_sparse ? '--sparse' : ''} \\
        -i salmon_index \\
        --gencode
    """
}
