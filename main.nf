#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { validateInput }           from './workflows/validate_input'
include { DISCOVERY }               from './workflows/discovery'
include { QUANTIFICATION }          from './workflows/quantification'
include { BUILD_DECOY_INDEX }       from './modules/local/build_decoy_index'
include { SALMON_QUANTIFY }         from './modules/local/salmon_quantify'
include { CHECK_STRANDEDNESS }      from './modules/local/check_strandedness'
include { CHECK_SALMON_LIBTYPE }    from './modules/local/check_salmon_libtype'
include { SPLIT_BY_NORMGROUP }      from './modules/local/split_normgroup'
include { PUBLISH_RUN_INFO }        from './modules/local/publish_run_info'
include { VERSION_DISCOVERY }       from './modules/local/version_tools'
include { VERSION_CODING }          from './modules/local/version_tools'
include { FIGURES }                 from './modules/local/figures'
include { RECOVERY_CHECK }          from './modules/local/recovery_check'

workflow {

    validateInput(
        params.input,
        params.genome,
        params.gtf,
        params.contrasts,
        params.design_formula
    )

    // ── Phase A — Discovery (all samples, condition-blind) ─────────────────
    ch_samples = validateInput.out.ch_samplesheet
    DISCOVERY(
        ch_samples,
        params.genome,
        params.gtf,
        params.outdir
    )

    // ── A2b: Strandedness guard — fail before Salmon if the empirical
    // orientation disagrees with the declared contract. Tiny task: on failure
    // the run aborts and any queued Salmon task is terminated.
    CHECK_STRANDEDNESS(DISCOVERY.out.strand_consensus)

    // ── B1: Build decoy-aware Salmon index (shared across all groups) ─────
    ch_transcripts_fa = channel.value(params.gencode_transcripts_fa ?: "")
    BUILD_DECOY_INDEX(
        DISCOVERY.out.frozen_gtf,
        params.genome,
        ch_transcripts_fa
    )
    ch_salmon_index = BUILD_DECOY_INDEX.out.index_dir.first()

    // ── B2: Salmon quantification — all samples, once ─────────────────────
    // Force a channel dependency so SALMON_QUANTIFY only starts after the
    // strandedness guard has passed (not merely "runs earlier").
    ch_samples_for_salmon = ch_samples
        .combine(CHECK_STRANDEDNESS.out.summary)
        .map { sample, cond, batch, fq1, fq2, bam, extra, _summary ->
            tuple(sample, cond, batch, fq1, fq2, bam, extra)
        }
    SALMON_QUANTIFY(ch_samples_for_salmon, ch_salmon_index)

    // ── B2b: Per-sample library-type guard (fails early on a wrong -l) ─────
    CHECK_SALMON_LIBTYPE(SALMON_QUANTIFY.out.quant_results)

    // Single shared value = full list of quant dirs; combine() keeps it as one
    // field instead of flattening the tuples into the group tuple
    ch_salmon_quants = CHECK_SALMON_LIBTYPE.out.checked
        .map { _sample, quant_dir -> quant_dir }
        .collect()
        .map { q -> [ q ] }

    // ── Build global contrasts TSV (used by all groups) ────────────────────
    ch_ct_lines = validateInput.out.ch_contrasts
        .map { name, num, den, batch ->
            "${name}\t${num}\t${den}\t${batch ?: ''}"
        }
    ch_contrasts_tsv = ch_ct_lines
        .collectFile(name: 'contrasts.tsv', newLine: true)

    // ── Split samplesheet by norm_group ────────────────────────────────────
    SPLIT_BY_NORMGROUP(
        validateInput.out.coldata_csv,
        validateInput.out.contrasts_csv
    )

    // Combine per-group coldata with shared Salmon quants + contrasts TSV
    ch_group_map = SPLIT_BY_NORMGROUP.out.per_group
        .flatten()
        .map { coldata_file ->
            def group_name = coldata_file.baseName
                .replaceAll('^coldata_', '')
                .replaceAll('\\.csv$', '')
            tuple(group_name, coldata_file)
        }
    ch_group_quants = ch_group_map
        .combine(ch_salmon_quants)
    ch_per_group = ch_group_quants
        .combine(ch_contrasts_tsv)

    // ── Phase B — Quantification & Statistics (per norm_group) ─────────────
    // DE runs on the full transcriptome (reference + novel) via the full-universe
    // tx2gene; origin (novel/annotated/reference_other) is added at reporting time.
    QUANTIFICATION(
        ch_per_group,
        DISCOVERY.out.tx2gene_detailed,
        validateInput.out.design_formula,
        DISCOVERY.out.gene_catalog,
        DISCOVERY.out.catalog_genes
    )

    // ── F4.6 / F4.7: run metadata + generated READMEs ──────────────────────
    // Samples index keyed by sample (norm_group comes from the covariates).
    ch_sample_rows = ch_samples.map { sample, cond, batch, _fq1, _fq2, _bam, extra ->
        "${sample}\t${extra.norm_group ?: 'NA'}\t${cond}\t${batch ?: 'NA'}"
    }
    ch_sample_index = ch_sample_rows
        .collectFile(name: 'sample_index.tsv', newLine: true,
                     seed: 'sample\tnorm_group\tcondition\tbatch\n')
        .first()

    // Version probes: one task per heavy image, once, so tool scripts stay cached.
    VERSION_DISCOVERY()
    VERSION_CODING()

    // F6.1/F6.2: reviewer validations (catalog + pooled expression). Cheap, once.
    ch_rec_tpm = QUANTIFICATION.out.expression_tpm
        .map { _group, f -> f }
        .collect()
    RECOVERY_CHECK(DISCOVERY.out.catalog_genes, ch_rec_tpm)

    // Every tool process writes a file named versions.yml, so the mixed channel
    // holds N homonymous paths and Nextflow refuses to stage them into one task
    // ("input file name collision"). Rename each file after its work-dir hash
    // (stable under -resume, so the consumer task stays cacheable) and collect.
    ch_versions = channel.empty()
        .mix(VERSION_DISCOVERY.out.versions,
             VERSION_CODING.out.versions,
             DISCOVERY.out.versions,
             QUANTIFICATION.out.versions,
             RECOVERY_CHECK.out.versions)
        .collectFile { f -> ["${f.getParent().baseName}.versions.yml", f.text] }
        .collect()

    PUBLISH_RUN_INFO(
        ch_sample_index,
        CHECK_STRANDEDNESS.out.summary,
        ch_versions,
        DISCOVERY.out.novelty,
        DISCOVERY.out.thresholds,
        DISCOVERY.out.gene_catalog,
        DISCOVERY.out.tier_sensitivity,
        QUANTIFICATION.out.qc_controls,
        QUANTIFICATION.out.universe,
        DISCOVERY.out.coding_consensus,
        DISCOVERY.out.recurrence,
        QUANTIFICATION.out.length_offset,
        channel.value(file("${projectDir}/docs/deliverables_contract.json"))
    )

    // F5.1 / F5.4: catalog figures (non-blocking, off with --run_figures false).
    // The expression-by-position panel reuses the per-group TPM tables plus the
    // sample index (condition + norm_group) from PUBLISH_RUN_INFO and the run's
    // comparisons CSV, whose denominators define the control conditions.
    if (params.run_figures) {
        ch_fig_tpm = QUANTIFICATION.out.expression_tpm
            .map { _group, f -> f }
            .collect()
        FIGURES(DISCOVERY.out.catalog_genes, DISCOVERY.out.tier_by_class,
                ch_fig_tpm, PUBLISH_RUN_INFO.out.sample_index,
                validateInput.out.contrasts_csv)
    }
}
