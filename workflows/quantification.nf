// Quantification & Statistics
// Receives frozen Discovery annotation + original FASTQs

include { BUILD_DECOY_INDEX }   from '../modules/local/build_decoy_index'
include { SALMON_QUANTIFY }     from '../modules/local/salmon_quantify'
include { TXIMPORT_AND_FILTER } from '../modules/local/tximport_filter'
include { DESEQ2_DE }           from '../modules/local/deseq2_de'
include { CIS_ASSOCIATIONS }    from '../modules/local/cis_associations'
include { FUNCTIONAL_ENRICHMENT } from '../modules/local/functional_enrichment'

workflow QUANTIFICATION {
    take:
    ch_samplesheet      // tuple: sample, condition, batch, fastq_1, fastq_2, bam, covariates
    frozen_gtf
    tx2gene_detailed
    ch_contrasts        // tuple: contrast_name, numerator, denominator, batch
    genome_fa
    reference_gtf       // original reference GTF (for cis-association gene lookup)
    design_formula
    outdir
    coldata_csv         // path to original samplesheet CSV (for R modules)

    main:
    // B1: Build decoy-aware Salmon index
    ch_salmon_index = BUILD_DECOY_INDEX(frozen_gtf, genome_fa, params.gencode_transcripts_fa)

    // B2: Salmon quantification per sample, collected for downstream
    ch_salmon_quants = SALMON_QUANTIFY(ch_samplesheet, ch_salmon_index)
        .map { sample, quant_dir -> tuple(sample, quant_dir) }
        .collect()

    // B3: tximport + independent filtering (condition-blind)
    TXIMPORT_AND_FILTER(
        ch_salmon_quants,
        tx2gene_detailed,
        coldata_csv
    )

    // Duplicate counts for DESeq2 + cis associations
    TXIMPORT_AND_FILTER.out.counts_rds.into { ch_counts_de; ch_counts_cis }

    // Write contrasts to a staged TSV for DESeq2
    ch_contrasts_tsv = ch_contrasts
        .map { name, num, den, batch ->
            "${name}\t${num}\t${den}\t${batch ?: ''}"
        }
        .collectFile(name: 'contrasts.tsv', newLine: true, seed: 'contrast_name\tnumerator\tdenominator\tbatch')

    // B4: DESeq2 differential expression
    DESEQ2_DE(
        ch_counts_de,
        coldata_csv,
        ch_contrasts_tsv,
        design_formula
    )

    // B5: Cis-regulatory candidate associations
    CIS_ASSOCIATIONS(
        DESEQ2_DE.out.contrasts_summary,
        ch_counts_cis,
        frozen_gtf,
        reference_gtf
    )

    // B6: Functional enrichment (ORA + GSEA)
    FUNCTIONAL_ENRICHMENT(
        CIS_ASSOCIATIONS.out.cis_pairs_sig,
        DESEQ2_DE.out.contrasts_summary
    )

    emit:
    deseq2_results   = DESEQ2_DE.out.all_results
    cis_results      = CIS_ASSOCIATIONS.out.cis_pairs_sig
    enrichment       = FUNCTIONAL_ENRICHMENT.out.ora_results
}
