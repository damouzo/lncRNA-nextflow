// Quantification & Statistics (per norm_group)
// Receives frozen Discovery annotation + a single norm_group's coldata + all Salmon quants
// Runs B3-B6 per norm_group to ensure independent size-factor normalisation.
//
// Ordering guarantee: every per-group artifact carries its norm_group as a
// tuple key (val(group), path(...)) and downstream stages pair channels with
// .join() on that key. Channels are NEVER matched positionally, because the
// emission order of a process output is completion order, which is not the
// submission order of the input channel.

include { TXIMPORT_AND_FILTER } from '../modules/local/tximport_filter'
include { DESEQ2_DE }           from '../modules/local/deseq2_de'
include { CIS_ASSOCIATIONS }    from '../modules/local/cis_associations'
include { FUNCTIONAL_ENRICHMENT } from '../modules/local/functional_enrichment'
include { REPORT_CONTRAST }      from '../modules/local/report_contrast'

workflow QUANTIFICATION {
    take:
    // ch_per_group: [group, coldata_csv, quant_tuples, contrasts_tsv_path]
    ch_per_group
    tx2gene_detailed
    design_formula
    gene_catalog
    outdir

    main:
    // B3: tximport + independent filtering (condition-blind, keyed by group)
    ch_txi_in = ch_per_group.map { group, coldata, quants, ct_tsv ->
        tuple(group, quants, coldata)
    }
    TXIMPORT_AND_FILTER(ch_txi_in, tx2gene_detailed)

    // B4: DESeq2 differential expression
    // Join the finished counts with each group's coldata + contrasts by key;
    // never positionally — completion order != submission order.
    ch_deseq_key = ch_per_group.map { group, coldata, quants, ct_tsv ->
        tuple(group, coldata, ct_tsv)
    }
    ch_deseq = TXIMPORT_AND_FILTER.out.counts_rds
        .join(ch_deseq_key)
        .map { group, counts, coldata, ct_tsv ->
            tuple(group, counts, coldata, ct_tsv)
        }

    DESEQ2_DE(ch_deseq, design_formula)

    // Summaries are already keyed by the process (val(group) in the tuple)
    ch_deseq_summary = DESEQ2_DE.out.contrasts_summary

    // B5: Cis-regulatory candidate associations (novel + annotated lncRNAs)
    ch_cis = TXIMPORT_AND_FILTER.out.counts_rds
        .join(ch_deseq_summary)
        .map { group, counts, summary ->
            tuple(group, counts, summary)
        }

    CIS_ASSOCIATIONS(ch_cis, gene_catalog)

    // B6: Functional enrichment (ORA + GSEA on protein-coding cis targets)
    ch_enrich = CIS_ASSOCIATIONS.out.cis_pairs_sig
        .join(ch_deseq_summary)
        .map { group, cisp, summary ->
            tuple(group, cisp, summary)
        }

    FUNCTIONAL_ENRICHMENT(ch_enrich, gene_catalog)

    // B7: Per-contrast reports (Quarto, self-contained HTML)
    // All phase-B outputs carry group as a tuple key, so the join is exact.
    // The per-group coldata and input contrasts file are passed so that:
    //   - every declared contrast that applies to THIS group gets a report even
    //     when it produced no results upstream, and
    //   - reports land flat (no per-group subdirectory).
    ch_ct_grp = ch_per_group.map { group, coldata, quants, ct_tsv ->
        tuple(group, coldata, ct_tsv)
    }
    ch_rpt = ch_deseq_summary
        .join(CIS_ASSOCIATIONS.out.cis_pairs)
        .join(CIS_ASSOCIATIONS.out.cis_pairs_sig)
        .join(FUNCTIONAL_ENRICHMENT.out.ora_results)
        .join(FUNCTIONAL_ENRICHMENT.out.gsea_results)
        .join(ch_ct_grp)
        .map { g, de, cisp, cissig, ora, gsea, coldata, ct_tsv ->
            tuple(g, de, cisp, cissig, ora, gsea, coldata, ct_tsv)
        }

    REPORT_CONTRAST(
        ch_rpt,
        Channel.value(file("${params.report_template ?: projectDir + '/assets/report_template.qmd'}")),
        gene_catalog
    )

    emit:
    deseq2_results   = DESEQ2_DE.out.all_results
    cis_results      = CIS_ASSOCIATIONS.out.cis_pairs_sig
    enrichment       = FUNCTIONAL_ENRICHMENT.out.ora_results
    reports          = REPORT_CONTRAST.out.reports
}