// Discovery: Assembly, Curation & Freezing
// Strictly condition-blind — no contrast-based filtering anywhere in this phase

include { BAM_QUICKCHECK }          from '../modules/local/bam_quickcheck'
include { BUILD_INFER_BED }         from '../modules/local/build_infer_bed'
include { INFER_STRANDEDNESS }      from '../modules/local/infer_strandedness'
include { STRINGTIE2_PER_SAMPLE }   from '../modules/local/stringtie_per_sample'
include { STRINGTIE2_MERGE }        from '../modules/local/stringtie_merge'
include { GFFCOMPARE }              from '../modules/local/gffcompare'
include { TRANSCRIPT_LENGTH_FILTER } from '../modules/local/length_filter'
include { CPAT_BUILD_MODEL }        from '../modules/local/cpat_build_model'
include { CPAT }                    from '../modules/local/cpat'
include { CPC2 }                    from '../modules/local/cpc2'
include { CODING_CONSENSUS }        from '../modules/local/coding_consensus'
include { CATALOG_OVERLAP }         from '../modules/local/catalog_overlap'
include { FEATURECOUNTS_EXPRESSION } from '../modules/local/featurecounts_expression'
include { EXPRESSION_RECURRENCE }   from '../modules/local/expression_recurrence'
include { ANNOTATION_FREEZE }       from '../modules/local/annotation_freeze'
include { BUILD_GENE_CATALOG }      from '../modules/local/build_gene_catalog'
include { CHECK_REFERENCE_COMPATIBILITY } from '../modules/local/check_reference_compatibility'
include { CONSERVATION_SCORE }      from '../modules/local/conservation_score'
include { SYNTENY_CHECK }           from '../modules/local/synteny_check'
include { ASSIGN_CATALOG_IDS }      from '../modules/local/assign_catalog_ids'
include { BUILD_CATALOG_MASTER }    from '../modules/local/build_catalog_master'
include { BUILD_F1_REGIONS }        from '../modules/local/build_f1_regions'
include { FLANK_COUNTS }            from '../modules/local/flank_counts'
include { FLANK_COVERAGE_RATIO }    from '../modules/local/flank_coverage_ratio'
include { ANTISENSE_COUNTS }        from '../modules/local/antisense_counts'
include { FLAG_ANTISENSE_LEAKAGE }  from '../modules/local/flag_antisense_leakage'
include { DERIVE_THRESHOLDS }       from '../modules/local/derive_thresholds'
include { COMPUTE_TIER_FINAL }      from '../modules/local/compute_tier_final'

workflow DISCOVERY {
    take:
    ch_samplesheet   // tuple: sample, condition, batch, fastq_1, fastq_2, bam, covariates
    genome_fa
    reference_gtf
    _outdir

    main:
    // Extract BAM paths for validation
    ch_bams = ch_samplesheet.map { sample, _cond, _batch, _fq1, _fq2, bam, _extra ->
        tuple(sample, bam)
    }

    // Validate BAM integrity
    BAM_QUICKCHECK(ch_bams)

    // External conservation/synteny resources must match params.genome.
    // Fail fast (header-only) before any bigWig/chain read — only runs when
    // at least one resource is configured. No-op otherwise (zero regression).
    // Resources are read from params inside the module (the reference dir is
    // mounted in the container, so no staged copy of the multi-GB files).
    if (params.conservation_bigwig || params.synteny_chain_file) {
        CHECK_REFERENCE_COMPATIBILITY(genome_fa)
    }

    // Duplicate validated BAMs: strandedness + StringTie2 + featureCounts
    ch_bams_strand   = BAM_QUICKCHECK.out.validated
    ch_bams_assembly = BAM_QUICKCHECK.out.validated
    ch_bams_fc       = BAM_QUICKCHECK.out.validated

    // Empirical strandedness detection. The reference BED12 is built once, not
    // per sample, and only the file paths are collected (CHECK_STRANDEDNESS reads
    // one summary table across all samples).
    BUILD_INFER_BED(reference_gtf)
    ch_strand_files = INFER_STRANDEDNESS(ch_bams_strand, BUILD_INFER_BED.out.bed12)
        .map { _sample, strandedness_txt -> strandedness_txt }
        .collect()

    // Per-sample StringTie2 assembly — collect GTFs for merge
    STRINGTIE2_PER_SAMPLE(ch_bams_assembly, reference_gtf)
    ch_gtf_list = STRINGTIE2_PER_SAMPLE.out.transcript_gtf
        .map { _sample, gtf -> gtf }
        .collect()

    // Merge StringTie2 assemblies
    STRINGTIE2_MERGE(ch_gtf_list, reference_gtf)

    // gffcompare classification vs reference
    GFFCOMPARE(STRINGTIE2_MERGE.out.merged_gtf, reference_gtf)

    // Length filter on annotated GTF
    TRANSCRIPT_LENGTH_FILTER(
        GFFCOMPARE.out.annotated_gtf,
        GFFCOMPARE.out.class_summary
    )

    // Duplicate length-filtered for CPAT, CPC2, featureCounts, and annotation_freeze
    ch_for_cpat   = TRANSCRIPT_LENGTH_FILTER.out.filtered_gtf
    ch_for_cpc2   = TRANSCRIPT_LENGTH_FILTER.out.filtered_gtf
    ch_for_fc     = TRANSCRIPT_LENGTH_FILTER.out.filtered_gtf
    ch_for_freeze = TRANSCRIPT_LENGTH_FILTER.out.filtered_gtf

    // Build CPAT models from reference GTF if not provided
    if (params.cpat_hexamer && params.cpat_logit_model) {
        ch_cpat_models = channel.value([
            file(params.cpat_hexamer, checkIfExists: true),
            file(params.cpat_logit_model, checkIfExists: true)
        ])
    } else {
        CPAT_BUILD_MODEL(genome_fa, reference_gtf)
        ch_cpat_models = CPAT_BUILD_MODEL.out.hexamer_table
            .combine(CPAT_BUILD_MODEL.out.logit_model)
            .map { hex, logit -> [hex, logit] }
    }

    CPAT(ch_for_cpat, genome_fa, ch_cpat_models)
    CPC2(ch_for_cpc2, genome_fa)

    CODING_CONSENSUS(CPAT.out.cpat_output, CPC2.out.cpc2_output, ch_for_cpc2)

    CATALOG_OVERLAP(
        CODING_CONSENSUS.out.consensus_table,
        CODING_CONSENSUS.out.consensus_gtf
    )

    // FeatureCounts: quantify all candidate transcripts against validated BAMs
    FEATURECOUNTS_EXPRESSION(ch_bams_fc, ch_for_fc)

    // Collect all per-sample featureCounts TSVs
    ch_fc_collected = FEATURECOUNTS_EXPRESSION.out.per_sample_counts
        .collect()

    // Expression recurrence filter (condition-blind) using FC counts
    EXPRESSION_RECURRENCE(
        CATALOG_OVERLAP.out.overlap_table,
        ch_fc_collected,
        TRANSCRIPT_LENGTH_FILTER.out.filtered_gtf
    )

    ANNOTATION_FREEZE(
        EXPRESSION_RECURRENCE.out.recurrence_table,
        ch_for_freeze,
        genome_fa,
        reference_gtf
    )

    // Conservation & synteny (reporting-only, optional). Gated so disabled
    // runs produce identical output; a header-only placeholder keeps
    // BUILD_GENE_CATALOG's join total (all-NA new columns) when skipped.
    if (params.conservation_bigwig) {
        CONSERVATION_SCORE(ANNOTATION_FREEZE.out.frozen_gtf)
        ch_conservation = CONSERVATION_SCORE.out.scores.first()
    } else {
        ch_conservation = channel.of("transcript_id\tmean_score\tmax_score\tpct_bases_conserved")
            .collectFile(name: 'conservation_placeholder.tsv', newLine: true).first()
    }

    if (params.synteny_chain_file && params.synteny_target_gtf) {
        SYNTENY_CHECK(ANNOTATION_FREEZE.out.frozen_gtf)
        ch_synteny = SYNTENY_CHECK.out.scores.first()
    } else {
        ch_synteny = channel.of("transcript_id\tsyntenic_locus\tsyntenic_target_gene_id")
            .collectFile(name: 'synteny_placeholder.tsv', newLine: true).first()
    }

    // Stable catalog identity: hierarchical IDs (gene + transcript) backed by a
    // persistent two-level registry outside results/. Runs AFTER freeezing so the
    // registry only ever sees transcripts that passed all discovery filters.
    // "cache false" inside the process is a deliberate -resume exception (it
    // writes outside the Nextflow work dir).
    ASSIGN_CATALOG_IDS(
        ANNOTATION_FREEZE.out.frozen_gtf,
        params.catalog_registry,        // val — may not exist on the first run
        params.genome_build,
        params.gtf,
        params.catalog_registry_project
    )

    // Gene catalog: one row per gene of the full analysis universe
    BUILD_GENE_CATALOG(
        ANNOTATION_FREEZE.out.frozen_gtf,
        reference_gtf,
        ch_conservation,
        ch_synteny,
        ASSIGN_CATALOG_IDS.out.catalog_ids,
        CATALOG_OVERLAP.out.overlap_table
    )

    // Novel catalog master (gene-level + transcript-level) with join assertions
    BUILD_CATALOG_MASTER(
        ASSIGN_CATALOG_IDS.out.catalog_ids,
        BUILD_GENE_CATALOG.out.gene_catalog,
        ASSIGN_CATALOG_IDS.out.catalog_tx,
        CATALOG_OVERLAP.out.overlap_table,
        EXPRESSION_RECURRENCE.out.recurrence_table,
        ANNOTATION_FREEZE.out.frozen_gtf
    )

    // ── F1: candidate/control regions and flags (Phase A, post-freeze) ─────
    // Terminal: nothing downstream depends on these outputs (design rule B).
    // Counts are per-sample (parallel); a single aggregation task merges the
    // matrices and computes metrics. Thresholds are a separate cheap process, so
    // the threshold rule can change with -resume without re-counting the BAMs.
    //
    // Channel cardinality: the single-file channels feeding the per-sample
    // counters are collapsed to value channels with .first(). A queue channel
    // paired with the 50-item BAM channel is zipped 1:1, which would silently
    // run the counter on one sample and derive thresholds from a partial cohort.
    ch_f1_coldata = ch_samplesheet
        .map { sample, _cond, _batch, _fq1, _fq2, bam, extra ->
            "${sample}\t${extra.norm_group ?: 'NA'}\t${bam.name}" }
        .collectFile(name: 'f1_coldata.tsv', newLine: true)
        .first()                                     // val — one cohort coldata

    BUILD_F1_REGIONS(
        ANNOTATION_FREEZE.out.frozen_gtf,
        reference_gtf,
        EXPRESSION_RECURRENCE.out.recurrence_table
    )
    // val — one region table and SAFs, reused by every sample's counters
    ch_f1_regions   = BUILD_F1_REGIONS.out.regions.first()
    ch_f1_flank_saf = BUILD_F1_REGIONS.out.flank_saf.first()
    ch_f1_anti_saf  = BUILD_F1_REGIONS.out.antisense_saf.first()

    // F1.3 (mono-exonic flank coverage) is exploratory and off by default: it did
    // not pass its pre-specified control separation, so it does not gate tier_final
    // (categorical mono-exonic rule replaces it). When disabled we feed an empty
    // flank_metrics table so the population stays unchanged and no flank work runs.
    if (params.run_flank_check) {
        FLANK_COUNTS(BAM_QUICKCHECK.out.validated, ch_f1_flank_saf)
        ch_flank_s2 = FLANK_COUNTS.out.s2.map { _sample, f -> f }.collect()
        ch_flank_s1 = FLANK_COUNTS.out.s1.map { _sample, f -> f }.collect()
        FLANK_COVERAGE_RATIO(
            ch_flank_s2,
            ch_flank_s1,
            ch_f1_coldata,
            ch_f1_regions
        )
        ch_flank_metrics = FLANK_COVERAGE_RATIO.out.flank_metrics
    } else {
        ch_flank_metrics = channel.of(
            "element_id\tset\tnorm_group\tflank_ratio_median\topp_strand_ratio_median\tn_samples_used\tflank_len_used"
        ).collectFile(name: 'flank_metrics_disabled.tsv', newLine: true).first()
        log.info "F1.3 flank check disabled (run_flank_check=false)"
    }

    ANTISENSE_COUNTS(BAM_QUICKCHECK.out.validated, ch_f1_anti_saf)
    ch_anti_s2 = ANTISENSE_COUNTS.out.s2.map { _sample, f -> f }.collect()
    ch_anti_s1 = ANTISENSE_COUNTS.out.s1.map { _sample, f -> f }.collect()
    FLAG_ANTISENSE_LEAKAGE(
        ch_anti_s2,
        ch_anti_s1,
        ch_f1_coldata,
        ch_f1_regions
    )

    DERIVE_THRESHOLDS(
        ch_flank_metrics,
        FLAG_ANTISENSE_LEAKAGE.out.antisense_metrics
    )
    COMPUTE_TIER_FINAL(
        ANNOTATION_FREEZE.out.frozen_gtf,
        reference_gtf,
        ch_flank_metrics,
        FLAG_ANTISENSE_LEAKAGE.out.antisense_metrics,
        DERIVE_THRESHOLDS.out.thresholds,
        DERIVE_THRESHOLDS.out.sensitivity
    )

    emit:
    frozen_gtf        = ANNOTATION_FREEZE.out.frozen_gtf
    tx2gene_detailed   = ANNOTATION_FREEZE.out.tx2gene_detailed
    gene_catalog       = BUILD_GENE_CATALOG.out.gene_catalog
    strand_consensus   = ch_strand_files
    catalog_ids        = ASSIGN_CATALOG_IDS.out.catalog_ids
    catalog_tx         = ASSIGN_CATALOG_IDS.out.catalog_tx
    catalog_gtf        = ASSIGN_CATALOG_IDS.out.catalog_gtf
    catalog_master     = BUILD_CATALOG_MASTER.out.master
    catalog_transcripts= BUILD_CATALOG_MASTER.out.transcripts
}
