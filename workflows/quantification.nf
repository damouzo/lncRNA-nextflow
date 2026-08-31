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
include { HETEROGENEITY_ANALYSIS } from '../modules/local/heterogeneity_analysis'
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
    // tx2gene_detailed and gene_catalog come from single-execution Discovery processes:
    // .first() turns them into reusable value channels so they pair with EVERY norm_group
    // instead of zipping positionally and silently limiting execution to one group.
    ch_tx2gene      = tx2gene_detailed.first()
    ch_gene_catalog = gene_catalog.first()

    // B3: tximport + independent filtering (condition-blind, keyed by group)
    ch_txi_in = ch_per_group.map { group, coldata, quants, ct_tsv ->
        tuple(group, quants, coldata)
    }
    TXIMPORT_AND_FILTER(ch_txi_in, ch_tx2gene)

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

    CIS_ASSOCIATIONS(ch_cis, ch_gene_catalog)

    // B6: Functional enrichment (ORA + GSEA on protein-coding cis targets)
    ch_enrich = CIS_ASSOCIATIONS.out.cis_pairs_sig
        .join(ch_deseq_summary)
        .map { group, cisp, summary ->
            tuple(group, cisp, summary)
        }

    FUNCTIONAL_ENRICHMENT(ch_enrich, ch_gene_catalog)

    // B5.5 (optional): Expression-heterogeneity diagnostic (Carrasco-Leon et al. 2021),
    // complementary to DESeq2. Opt-in per norm_group via params.heterogeneity_analysis;
    // groups not listed there default to disabled (fail-safe). The input channel is
    // filtered by an inner .join() against the enabled-groups channel, so the process is
    // never scheduled at all for a disabled group (not submitted-then-skipped).
    ch_het_enabled_groups = ch_per_group
        .map { group, coldata, quants, ct_tsv -> group }
        .filter { group -> (params.heterogeneity_analysis?.get(group) ?: false) as boolean }
        .map { group -> tuple(group, true) }

    ch_het_key = ch_per_group.map { group, coldata, quants, ct_tsv ->
        tuple(group, coldata, ct_tsv)
    }
    ch_het_in = TXIMPORT_AND_FILTER.out.counts_rds
        .join(TXIMPORT_AND_FILTER.out.tpm_rds)
        .join(ch_het_key)
        .join(ch_het_enabled_groups)
        .map { group, counts, tpm, coldata, ct_tsv, flag ->
            tuple(group, counts, tpm, coldata, ct_tsv)
        }

    HETEROGENEITY_ANALYSIS(ch_het_in, ch_gene_catalog)

    // Disabled groups still need one channel entry each so the REPORT_CONTRAST join
    // below never drops them — a header-only placeholder (0 data rows) is broadcast to
    // every disabled group. REPORT_CONTRAST treats a 0-row heterogeneity table the same
    // way it already treats other empty inputs: the new report section is simply skipped,
    // so disabled-group reports render identically to the pre-feature pipeline.
    ch_het_disabled_groups = ch_per_group
        .map { group, coldata, quants, ct_tsv -> group }
        .filter { group -> !((params.heterogeneity_analysis?.get(group) ?: false) as boolean) }

    ch_het_gene_placeholder = Channel
        .of("contrast\tgene_id\tn_case\tn_control\tn_up\tn_down\tn_no_change\tconcordant_frac\tdiscordant_frac\tdirection\taltered_call")
        .collectFile(name: 'heterogeneity_placeholder.tsv', newLine: true)
        .first()
    ch_het_cv_placeholder = Channel
        .of("contrast\tgene_id\torigin_group\tcv\tn_lncrna_genes\tn_coding_genes\tmean_cv_lncrna\tmean_cv_coding\tmedian_cv_lncrna\tmedian_cv_coding\tt_statistic\tp_value")
        .collectFile(name: 'heterogeneity_cv_placeholder.tsv', newLine: true)
        .first()

    ch_het_summary_all = HETEROGENEITY_ANALYSIS.out.heterogeneity_summary
        .mix(ch_het_disabled_groups.combine(ch_het_gene_placeholder))
    ch_het_cv_summary_all = HETEROGENEITY_ANALYSIS.out.heterogeneity_cv_summary
        .mix(ch_het_disabled_groups.combine(ch_het_cv_placeholder))

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
        .join(ch_het_summary_all)
        .join(ch_het_cv_summary_all)
        .join(ch_ct_grp)
        .join(TXIMPORT_AND_FILTER.out.tpm_rds)
        .map { g, de, cisp, cissig, ora, gsea, het, hetcv, coldata, ct_tsv, tpm ->
            tuple(g, de, cisp, cissig, ora, gsea, het, hetcv, coldata, ct_tsv, tpm)
        }

    REPORT_CONTRAST(
        ch_rpt,
        Channel.value(file("${params.report_template ?: projectDir + '/assets/report_template.qmd'}")),
        ch_gene_catalog
    )

    emit:
    deseq2_results   = DESEQ2_DE.out.all_results
    cis_results      = CIS_ASSOCIATIONS.out.cis_pairs_sig
    enrichment       = FUNCTIONAL_ENRICHMENT.out.ora_results
    heterogeneity    = HETEROGENEITY_ANALYSIS.out.heterogeneity_summary
    reports          = REPORT_CONTRAST.out.reports
}