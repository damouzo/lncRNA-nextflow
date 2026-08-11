// Quantification & Statistics (per norm_group)
// Receives frozen Discovery annotation + a single norm_group's coldata + all Salmon quants
// Runs B3-B6 per norm_group to ensure independent size-factor normalisation

include { TXIMPORT_AND_FILTER } from '../modules/local/tximport_filter'
include { DESEQ2_DE }           from '../modules/local/deseq2_de'
include { CIS_ASSOCIATIONS }    from '../modules/local/cis_associations'
include { FUNCTIONAL_ENRICHMENT } from '../modules/local/functional_enrichment'

workflow QUANTIFICATION {
    take:
    // ch_per_group: [group, coldata_csv, quant_tuples, contrasts_tsv_path]
    ch_per_group
    tx2gene_detailed
    design_formula
    frozen_gtf
    reference_gtf
    outdir

    main:
    ch_per_group
        .multiMap { group, coldata, quants, ct_tsv ->
            txi_quants:    tuple(quants)
            txi_coldata:   coldata
            txi_group:     group
            de_coldata:    coldata
            de_ctsv:       ct_tsv
            de_group:      group
            cis_group:     group
            enrich_group:  group
        }
        .set { ch_split }

    // B3: tximport + independent filtering (condition-blind, per-group coldata)
    TXIMPORT_AND_FILTER(
        ch_split.txi_quants,
        tx2gene_detailed,
        ch_split.txi_coldata,
        ch_split.txi_group
    )

    // Duplicate counts for DESeq2 + cis associations
    TXIMPORT_AND_FILTER.out.counts_rds.into { ch_counts_de; ch_counts_cis }

    // B4: DESeq2 differential expression
    DESEQ2_DE(
        ch_counts_de,
        ch_split.de_coldata,
        ch_split.de_ctsv,
        design_formula,
        ch_split.de_group
    )

    // Duplicate contrast summary for cis + enrichment
    DESEQ2_DE.out.contrasts_summary.into { ch_de_summary_cis; ch_de_summary_enrich }

    // B5: Cis-regulatory candidate associations
    CIS_ASSOCIATIONS(
        ch_de_summary_cis,
        ch_counts_cis,
        frozen_gtf,
        reference_gtf,
        ch_split.cis_group
    )

    // B6: Functional enrichment (ORA + GSEA)
    FUNCTIONAL_ENRICHMENT(
        CIS_ASSOCIATIONS.out.cis_pairs_sig,
        ch_de_summary_enrich,
        ch_split.enrich_group
    )

    emit:
    deseq2_results   = DESEQ2_DE.out.all_results
    cis_results      = CIS_ASSOCIATIONS.out.cis_pairs_sig
    enrichment       = FUNCTIONAL_ENRICHMENT.out.ora_results
}