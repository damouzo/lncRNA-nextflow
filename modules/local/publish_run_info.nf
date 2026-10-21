// F4.6 / F4.7: publish run metadata and the contract-derived README.
// Writes the audit set (params_used, versions, sample_index, run_manifest,
// external_catalogs) and the run QA report, plus the two generated READMEs.
// The audit files go to <outdir>/run_info/ (pattern-scoped in conf/modules.config);
// versions/sample_index/qa_report/README go to <outdir>/deliverables[/run_info]/.

process PUBLISH_RUN_INFO {

    label 'process_single'

    input:
    path samples              // sample,norm_group,condition,batch TSV
    path strandedness         // CHECK_STRANDEDNESS strandedness_summary.tsv
    path versions_files       // collected versions.yml from every tool process
    path novelty_status       // CHECK_NOVELTY novelty_status.tsv
    path thresholds           // DERIVE_THRESHOLDS thresholds_used.tsv
    path gene_catalog         // BUILD_GENE_CATALOG gene_catalog.tsv
    path sensitivity          // COMPUTE_TIER_FINAL tier_final_sensitivity.tsv
    path contract             // docs/deliverables_contract.json

    output:
    path "params_used.tsv"
    path "versions.tsv"
    path "sample_index.tsv"
    path "run_manifest.tsv"
    path "external_catalogs.tsv"
    path "qa_report.tsv"
    path "__deliverables_README.md"
    path "__run_info_README.md"

    script:
    // Names come from docs/external_catalogs.json (single source); paths from params.
    def sources = new groovy.json.JsonSlurper()
        .parse(new File("${projectDir}/docs/external_catalogs.json"))['sources']
    def cats = sources.collect { s ->
        def path = params["${s.name}_bed"]
        path ? "--catalog ${s.name}=\"${path}\"".toString() : null
    }.findAll { it != null }
    def cat_args = cats.join(' ')
    def ver_args = versions_files.collect { "--versions-yml \"${it}\"" }.join(' ')
    """
    build_run_info.py \\
        --outdir . \\
        --samples "${samples}" \\
        --strandedness "${strandedness}" \\
        ${ver_args} \\
        --novelty-status "${novelty_status}" \\
        --thresholds "${thresholds}" \\
        --gene-catalog "${gene_catalog}" \\
        --sensitivity "${sensitivity}" \\
        --pipeline-version "${workflow.manifest.version}" \\
        --publish-level "${params.publish_level}" \\
        --outdir-param "${params.outdir}" \\
        --species "${params.species}" \\
        --genome-build "${params.genome_build ?: ''}" \\
        --library-strandedness "${params.library_strandedness}" \\
        --design-formula "${params.design_formula}" \\
        --tier1-classes '${params.tier1_classes.join(",")}' \\
        --min-transcript-length "${params.min_transcript_length}" \\
        --fdr "${params.fdr_threshold}" \\
        --lfc "${params.lfc_threshold}" \\
        ${cat_args}

    build_deliverables_readme.py \\
        --contract "${contract}" \\
        --outdir . \\
        --deliverables-readme __deliverables_README.md \\
        --run-info-readme __run_info_README.md \\
        --qa-report qa_report.tsv \\
        --publish-level "${params.publish_level}" \\
        --fdr "${params.fdr_threshold}" \\
        --lfc "${params.lfc_threshold}"
    """
}
