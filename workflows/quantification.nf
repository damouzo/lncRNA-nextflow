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
include { SPLIT_CONTRASTS }      from '../modules/local/split_contrasts'
include { DESIGN_CHECK }         from '../modules/local/design_check'
include { MAKE_XLSX }            from '../modules/local/make_xlsx'
include { BUILD_CONTRAST_TABLES } from '../modules/local/build_contrast_tables'
include { SPLIT_ENRICHMENT }     from '../modules/local/split_enrichment'
include { BUILD_CONTRASTS_OVERVIEW } from '../modules/local/build_contrasts_overview'
include { CHECK_DELIVERABLES }   from '../modules/local/check_deliverables'
include { QC_CONTROLS }          from '../modules/local/qc_controls'

workflow QUANTIFICATION {
    take:
    // ch_per_group: [group, coldata_csv, quant_tuples, contrasts_tsv_path]
    ch_per_group
    tx2gene_detailed
    design_formula
    gene_catalog
    catalog_genes       // deliverables/catalog catalog_genes.tsv (annotation for de_full)
    novelty_status      // CHECK_NOVELTY novelty_status.tsv (run-wide, single value)
    strandedness        // CHECK_STRANDEDNESS strandedness_summary.tsv (run-wide)

    main:
    // tx2gene_detailed and gene_catalog come from single-execution Discovery processes:
    // .first() turns them into reusable value channels so they pair with EVERY norm_group
    // instead of zipping positionally and silently limiting execution to one group.
    ch_tx2gene      = tx2gene_detailed.first()
    ch_gene_catalog = gene_catalog.first()
    // catalog_genes.tsv is the annotated superset of gene_catalog.tsv (adds
    // lnc_position_class + gene_group); used by tximport TSVs, cis and the report.
    ch_catalog_genes = catalog_genes.first()
    ch_novelty      = novelty_status.first()
    ch_strandedness = strandedness.first()

    // (group, coldata, contrasts) — reused by DESIGN_CHECK and SPLIT_CONTRASTS.
    ch_ct_in = ch_per_group.map { group, coldata, _quants, ct_tsv ->
        tuple(group, coldata, ct_tsv)
    }

    // F3.2: per-contrast design rank check (metadata only, cheap, before DESeq2)
    DESIGN_CHECK(ch_ct_in, design_formula)

    // B3: tximport + independent filtering (condition-blind, keyed by group)
    ch_txi_in = ch_per_group.map { group, coldata, quants, _ct_tsv ->
        tuple(group, quants, coldata)
    }
    TXIMPORT_AND_FILTER(ch_txi_in, ch_tx2gene, ch_catalog_genes)

    // B4: DESeq2 differential expression
    // Join the finished counts with the ALIGNED coldata (counts column order) plus
    // contrasts by key. H0: never use SPLIT_BY_NORMGROUP's coldata here — its rows
    // follow the samplesheet, not the quant-dir order of the counts matrix.
    ch_deseq_key = ch_per_group.map { group, _coldata, _quants, ct_tsv ->
        tuple(group, ct_tsv)
    }
    ch_deseq = TXIMPORT_AND_FILTER.out.counts_rds
        .join(TXIMPORT_AND_FILTER.out.dds_rds)
        .join(TXIMPORT_AND_FILTER.out.aligned_coldata)
        .join(ch_deseq_key)
        .map { group, counts, dds, coldata, ct_tsv ->
            tuple(group, counts, dds, coldata, ct_tsv)
        }

    DESEQ2_DE(ch_deseq, design_formula)

    // Slice D: optional, non-blocking xlsx export of each contrast's DE table.
    ch_de_ct = DESEQ2_DE.out.de_per_contrast
        .flatMap { group, files ->
            def list = files instanceof List ? files : [files]
            list.collect { f ->
                tuple(group, f.name.replaceFirst(/_DE\.tsv$/, ''), f)
            }
        }

    // DDX41/DHX34 positive-control report over all DE tables (diagnostic only,
    // never fails the run). The manifest carries the
    // (group, contrast) of each table because a staged `path` list is flattened.
    // The process is launched at the very end (after BUILD_CONTRASTS_OVERVIEW) so
    // the guard runs after the deliverables are computed and published.
    ch_contrasts_ref = ch_per_group.map { _g, _c, _quants, ct_tsv -> ct_tsv }.first()
    ch_qc_de = ch_de_ct.map { _g, _ct, de_tsv -> de_tsv }.collect()
    ch_qc_manifest = ch_de_ct
        .map { g, ct, de_tsv -> "${g}\t${ct}\t${de_tsv.name}" }
        .collectFile(name: 'qc_de_manifest.tsv', newLine: true,
                     seed: 'group\tcontrast\tfile\n')

    if (params.export_xlsx) {
        MAKE_XLSX(ch_de_ct)
    }

    // Summaries are already keyed by the process (val(group) in the tuple)
    ch_deseq_summary = DESEQ2_DE.out.contrasts_summary

    // B5: Cis-regulatory candidate associations (novel + annotated lncRNAs)
    // H0: use the counts-aligned coldata, not SPLIT_BY_NORMGROUP's.
    ch_cis = TXIMPORT_AND_FILTER.out.counts_rds
        .join(ch_deseq_summary)
        .join(TXIMPORT_AND_FILTER.out.aligned_coldata)
        .map { group, counts, summary, coldata ->
            tuple(group, counts, summary, coldata)
        }

    CIS_ASSOCIATIONS(ch_cis, ch_catalog_genes)

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
        .map { group, _coldata, _quants, _ct_tsv -> group }
        .filter { group -> (params.heterogeneity_analysis?.get(group) ?: false) as boolean }
        .map { group -> tuple(group, true) }

    ch_het_key = ch_per_group.map { group, _coldata, _quants, ct_tsv ->
        tuple(group, ct_tsv)
    }
    ch_het_in = TXIMPORT_AND_FILTER.out.counts_rds
        .join(TXIMPORT_AND_FILTER.out.tpm_rds)
        .join(TXIMPORT_AND_FILTER.out.aligned_coldata)
        .join(ch_het_key)
        .join(ch_het_enabled_groups)
        .map { group, counts, tpm, coldata, ct_tsv, _flag ->
            tuple(group, counts, tpm, coldata, ct_tsv)
        }

    HETEROGENEITY_ANALYSIS(ch_het_in, ch_gene_catalog)

    // Disabled groups still need one channel entry each so the REPORT_CONTRAST join
    // below never drops them — a header-only placeholder (0 data rows) is broadcast to
    // every disabled group. REPORT_CONTRAST treats a 0-row heterogeneity table the same
    // way it already treats other empty inputs: the new report section is simply skipped,
    // so disabled-group reports render identically to the pre-feature pipeline.
    ch_het_disabled_groups = ch_per_group
        .map { group, _coldata, _quants, _ct_tsv -> group }
        .filter { group -> !((params.heterogeneity_analysis?.get(group) ?: false) as boolean) }

    ch_het_gene_placeholder = channel
        .of("contrast\tgene_id\tn_case\tn_control\tn_up\tn_down\tn_no_change\tconcordant_frac\tdiscordant_frac\tdirection\taltered_call")
        .collectFile(name: 'heterogeneity_placeholder.tsv', newLine: true)
        .first()
    ch_het_cv_placeholder = channel
        .of("contrast\tgene_id\torigin_group\tcv\tn_lncrna_genes\tn_coding_genes\tmean_cv_lncrna\tmean_cv_coding\tmedian_cv_lncrna\tmedian_cv_coding\tt_statistic\tp_value")
        .collectFile(name: 'heterogeneity_cv_placeholder.tsv', newLine: true)
        .first()

    ch_het_summary_all = HETEROGENEITY_ANALYSIS.out.heterogeneity_summary
        .mix(ch_het_disabled_groups.combine(ch_het_gene_placeholder))
    ch_het_cv_summary_all = HETEROGENEITY_ANALYSIS.out.heterogeneity_cv_summary
        .mix(ch_het_disabled_groups.combine(ch_het_cv_placeholder))

    // B7: Per-contrast reports (Quarto, self-contained HTML).
    // SPLIT_CONTRASTS writes one 1-row contrasts file per (group, contrast) that
    // applies to the group (ch_ct_in is defined at the top of main); each
    // REPORT_CONTRAST task then renders exactly that contrast, so publishDir can
    // land it in contrasts/<grupo>__<contraste>/ from the task inputs.
    SPLIT_CONTRASTS(ch_ct_in)

    ch_ct_keyed = SPLIT_CONTRASTS.out.per_contrast
        .flatMap { group, files ->
            def list = files instanceof List ? files : [files]
            def prefix = group + '__'
            list.collect { f ->
                def base = f.name
                if (!base.startsWith(prefix)) {
                    error("SPLIT_CONTRASTS: unexpected file ${base} for group ${group}")
                }
                tuple(group, base.substring(prefix.length()).replaceFirst(/\.tsv$/, ''), f)
            }
        }

    // All phase-B outputs carry group as a tuple key, so every join is exact.
    ch_rpt_base = ch_deseq_summary
        .join(CIS_ASSOCIATIONS.out.cis_pairs)
        .join(CIS_ASSOCIATIONS.out.cis_pairs_sig)
        .join(FUNCTIONAL_ENRICHMENT.out.ora_results)
        .join(FUNCTIONAL_ENRICHMENT.out.gsea_results)
        .join(ch_het_summary_all)
        .join(ch_het_cv_summary_all)
        .join(TXIMPORT_AND_FILTER.out.aligned_coldata)
        .join(TXIMPORT_AND_FILTER.out.tpm_rds)
        .map { g, de, cisp, cissig, ora, gsea, het, hetcv, coldata, tpm ->
            tuple(g, de, cisp, cissig, ora, gsea, het, hetcv, coldata, tpm)
        }

    // combine(by: 0), NOT join: ch_rpt_base is 1 item/group while ch_ct_keyed is
    // N items/group. JoinOp emits a single pair per key and silently drops the rest
    // under remainder=false; CombineOp with a pivot emits the full per-key product.
    // Output field order is (key, left-rest, right-rest); re-map explicitly to be
    // independent of which side emits first.
    ch_rpt = ch_rpt_base
        .combine(ch_ct_keyed, by: 0)
        .map { g, de, cisp, cissig, ora, gsea, het, hetcv, coldata, tpm, ct, ct_tsv ->
            tuple(g, ct, de, cisp, cissig, ora, gsea, het, hetcv, coldata, ct_tsv, tpm)
        }

    REPORT_CONTRAST(
        ch_rpt,
        channel.value(file("${params.report_template ?: projectDir + '/assets/report_template.qmd'}")),
        ch_catalog_genes,
        ch_novelty,
        ch_strandedness
    )

    // F4.2: per-contrast master tables — last node of the graph.
    // ch_de_ct is N items/group and the universe is 1/group, so this must be a
    // per-key product (combine), not a join (which would keep only the first).
    // The universe comes from TXIMPORT's pre-filter universe so nothing is lost.
    ch_bt_in = ch_de_ct
        .combine(TXIMPORT_AND_FILTER.out.universe_tsv.map { g, u -> tuple(g, u) }, by: 0)
        .map { g, ct, de_tsv, universe_tsv ->
            tuple(g, ct, de_tsv, universe_tsv)
        }
    BUILD_CONTRAST_TABLES(ch_bt_in, catalog_genes.first())

    // Split the group-level enrichment table (one file per group, rows keyed by the
    // `contrast` column) into one file per (group, contrast) so it is delivered next
    // to the DE tables it belongs to. Purely a repartition: no rows are duplicated
    // or dropped, and a header is always written.
    ch_split_enrich_in = FUNCTIONAL_ENRICHMENT.out.ora_results
        .join(FUNCTIONAL_ENRICHMENT.out.gsea_results)
        .combine(ch_ct_keyed.map { g, ct, _tsv -> tuple(g, ct) }, by: 0)
        .map { g, ora, gsea, ct ->
            tuple(g, ct, ora, gsea)
        }
    SPLIT_ENRICHMENT(ch_split_enrich_in)

    // One row per contrast: concatenation of every summary.tsv plus the count of
    // novel Tier-1 significant genes (join of de_full_significant with the catalog).
    // collectFile renames each file after its contrast id so the single downstream
    // task can stage them without the homonymous-name collision.
    ch_overview_summaries = BUILD_CONTRAST_TABLES.out.overview_input
        .collectFile { id, summary, _sig -> ["${id}__summary.tsv", summary.text] }
        .collect()
    ch_overview_significant = BUILD_CONTRAST_TABLES.out.overview_input
        .collectFile { id, _summary, sig -> ["${id}__significant.tsv", sig.text] }
        .collect()
    BUILD_CONTRASTS_OVERVIEW(ch_overview_summaries, ch_overview_significant,
                             catalog_genes.first())

    // Positive-control report, launched last: it depends on the overview so it
    // runs after every deliverables file (counts, DE tables, reports, overview)
    // has been computed and published. Diagnostic only.
    QC_CONTROLS(ch_qc_de, ch_qc_manifest, ch_contrasts_ref,
                BUILD_CONTRASTS_OVERVIEW.out.overview)

    // Final guard: fail the run if the per-contrast stage did not receive every
    // (group, contrast) that SPLIT_CONTRASTS declared for a group. The ids are
    // compared on the channel values, not on the published tree: publishDir copies
    // are asynchronous (FilePublisher submits them to a thread pool after outputs
    // are bound), so the delivered tree cannot be read reliably within the run.
    ch_expected_contrast_ids = ch_ct_keyed
        .map { g, ct, _tsv -> "${g}__${ct}" }
        .collect()
    ch_delivered_contrast_ids = ch_rpt
        .map { g, ct, _de, _cisp, _cissig, _ora, _gsea, _het, _hetcv, _coldata, _ct_tsv, _tpm ->
            "${g}__${ct}"
        }
        .collect()
    CHECK_DELIVERABLES(
        ch_expected_contrast_ids,
        ch_delivered_contrast_ids,
        ch_contrasts_ref,
        BUILD_CONTRASTS_OVERVIEW.out.overview,
        channel.value(file("${projectDir}/docs/deliverables_contract.json"))
    )

    ch_versions = channel.empty()
        .mix(TXIMPORT_AND_FILTER.out.versions,
             DESEQ2_DE.out.versions,
             CIS_ASSOCIATIONS.out.versions,
             FUNCTIONAL_ENRICHMENT.out.versions,
             HETEROGENEITY_ANALYSIS.out.versions,
             REPORT_CONTRAST.out.versions)
        .mix(BUILD_CONTRAST_TABLES.out.versions)
        .mix(DESIGN_CHECK.out.versions)
        .mix(SPLIT_ENRICHMENT.out.versions)
        .mix(BUILD_CONTRASTS_OVERVIEW.out.versions)
        .mix(CHECK_DELIVERABLES.out.versions)
        .mix(QC_CONTROLS.out.versions)

    emit:
    versions         = ch_versions
    expression_tpm   = TXIMPORT_AND_FILTER.out.tpm_tsv   // (group, <group>_tpm.tsv) for FIGURES
    deseq2_results   = DESEQ2_DE.out.all_results
    cis_results      = CIS_ASSOCIATIONS.out.cis_pairs_sig
    enrichment       = FUNCTIONAL_ENRICHMENT.out.ora_results
    heterogeneity    = HETEROGENEITY_ANALYSIS.out.heterogeneity_summary
    reports          = REPORT_CONTRAST.out.reports
    qc_controls      = QC_CONTROLS.out.qc_controls
    universe         = TXIMPORT_AND_FILTER.out.universe_tsv.map { _g, f -> f }.collect()
    length_offset    = DESEQ2_DE.out.length_offset.map { _g, f -> f }.collect()
}