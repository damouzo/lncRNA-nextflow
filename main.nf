#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { validateInput }           from './workflows/validate_input'
include { DISCOVERY }               from './workflows/discovery'
include { QUANTIFICATION }          from './workflows/quantification'

workflow {

    validateInput(
        params.input,
        params.genome,
        params.gtf,
        params.contrasts,
        params.design_formula
    )

    // Duplicate samplesheet channel — one copy for Discovery, one for Quantification
    validateInput.out.ch_samplesheet.into { ch_for_discovery; ch_for_quant }

    DISCOVERY(
        ch_for_discovery,
        params.genome,
        params.gtf,
        params.outdir
    )

    QUANTIFICATION(
        ch_for_quant,
        DISCOVERY.out.frozen_gtf,
        DISCOVERY.out.tx2gene_detailed,
        validateInput.out.ch_contrasts,
        params.genome,
        params.gtf,
        validateInput.out.design_formula,
        params.outdir,
        validateInput.out.coldata_csv
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
