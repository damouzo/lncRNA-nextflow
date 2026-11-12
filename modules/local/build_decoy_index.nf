// Build decoy-aware Salmon index combining the quantification universe + reference
// transcriptome. The universe carries the catalog transcripts AND the reference-
// derived novel isoforms (isoform_classes), so transcript-level EM can split reads
// between the reference isoforms and the novel isoform; BUILD_TX2GENE then sums
// those reads back under the reference parent gene.
// If gencode_transcripts_fa is not provided, derives the reference transcriptome
// from the GTF annotation + genome FASTA via gffread.

process BUILD_DECOY_INDEX {

    input:
    path quant_universe_gtf    // ANNOTATION_FREEZE quant_universe.gtf (catalog + isoforms)
    path genome_fa
    val  gencode_transcripts_fa

    output:
    path "salmon_index", emit: index_dir

    script:
    """
    mkdir -p salmon_index

    # Extract the novel (catalog + isoform) sequences from the quantification universe
    gffread -w lncrna_transcripts.fa -g "${genome_fa}" "${quant_universe_gtf}"

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
