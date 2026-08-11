#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { validateInput }           from './workflows/validate_input'
include { DISCOVERY }               from './workflows/discovery'
include { QUANTIFICATION }          from './workflows/quantification'
include { BUILD_DECOY_INDEX }       from './modules/local/build_decoy_index'
include { SALMON_QUANTIFY }         from './modules/local/salmon_quantify'
include { SPLIT_BY_NORMGROUP }      from './modules/local/split_normgroup'

workflow {

    validateInput(
        params.input,
        params.genome,
        params.gtf,
        params.contrasts,
        params.design_formula
    )

    // Duplicate samplesheet: Discovery + Salmon (both use all samples)
    validateInput.out.ch_samplesheet.into { ch_for_discovery; ch_for_salmon }

    // ── Phase A — Discovery (all samples, condition-blind) ─────────────────
    DISCOVERY(
        ch_for_discovery,
        params.genome,
        params.gtf,
        params.outdir
    )

    // ── B1: Build decoy-aware Salmon index (shared across all groups) ─────
    ch_salmon_index = BUILD_DECOY_INDEX(
        DISCOVERY.out.frozen_gtf,
        params.genome,
        params.gencode_transcripts_fa
    )

    // ── B2: Salmon quantification — all samples, once ─────────────────────
    ch_salmon_quants = SALMON_QUANTIFY(ch_for_salmon, ch_salmon_index)
        .map { sample, quant_dir -> tuple(sample, quant_dir) }
        .collect()

    // ── Build global contrasts TSV (used by all groups) ────────────────────
    ch_contrasts_tsv = validateInput.out.ch_contrasts
        .map { name, num, den, batch ->
            "${name}\t${num}\t${den}\t${batch ?: ''}"
        }
        .collectFile(name: 'contrasts.tsv', newLine: true,
                     seed: 'contrast_name\tnumerator\tdenominator\tbatch')

    // ── Split samplesheet by norm_group ────────────────────────────────────
    SPLIT_BY_NORMGROUP(
        validateInput.out.coldata_csv,
        validateInput.out.contrasts_csv
    )

    // Combine per-group coldata with shared Salmon quants + contrasts TSV
    ch_per_group = SPLIT_BY_NORMGROUP.out.per_group
        .map { coldata_file ->
            def group_name = coldata_file.baseName
                .replaceAll('^coldata_', '')
                .replaceAll('\\.csv$', '')
            tuple(group_name, coldata_file)
        }
        .combine(ch_salmon_quants)
        .combine(ch_contrasts_tsv)

    // ── Phase B — Quantification & Statistics (per norm_group) ─────────────
    QUANTIFICATION(
        ch_per_group,
        DISCOVERY.out.tx2gene_detailed,
        validateInput.out.design_formula,
        DISCOVERY.out.frozen_gtf,
        params.gtf,
        params.outdir
    )
}

workflow.onComplete {
    log.info """
    ─────────────────────────────────────────
    lncRNA-discovery pipeline finished
    Results: ${params.outdir}
    ─────────────────────────────────────────
    """.stripIndent()
}