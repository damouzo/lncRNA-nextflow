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

    // ── Phase A — Discovery (all samples, condition-blind) ─────────────────
    ch_samples = validateInput.out.ch_samplesheet
    ch_samples.count().view { "N_SAMPLES_REALES: $it" }
    DISCOVERY(
        ch_samples,
        params.genome,
        params.gtf,
        params.outdir
    )

    // ── B1: Build decoy-aware Salmon index (shared across all groups) ─────
    ch_transcripts_fa = Channel.value(params.gencode_transcripts_fa ?: "")
    BUILD_DECOY_INDEX(
        DISCOVERY.out.frozen_gtf,
        params.genome,
        ch_transcripts_fa
    )
    ch_salmon_index = BUILD_DECOY_INDEX.out.index_dir.first()

    // ── B2: Salmon quantification — all samples, once ─────────────────────
    SALMON_QUANTIFY(ch_samples, ch_salmon_index)
    // Single shared value = full list of quant dirs; combine() keeps it as one
    // field instead of flattening the tuples into the group tuple
    ch_salmon_quants = SALMON_QUANTIFY.out.quant_results
        .map { sample, quant_dir -> quant_dir }
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
    // TEMPORARY diagnostic: confirm all norm_groups reach Phase B (remove after validation)
    ch_group_map.map { g, f -> g }.collect().view { "NORM_GROUPS_FOR_PHASE_B: $it" }
    ch_group_quants = ch_group_map
        .combine(ch_salmon_quants)
    ch_per_group = ch_group_quants
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
