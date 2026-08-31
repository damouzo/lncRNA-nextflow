// Discovery: Assembly, Curation & Freezing
// Strictly condition-blind — no contrast-based filtering anywhere in this phase

include { BAM_QUICKCHECK }          from '../modules/local/bam_quickcheck'
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

workflow DISCOVERY {
    take:
    ch_samplesheet   // tuple: sample, condition, batch, fastq_1, fastq_2, bam, covariates
    genome_fa
    reference_gtf
    outdir

    main:
    // Extract BAM paths for validation
    ch_bams = ch_samplesheet.map { sample, cond, batch, fq1, fq2, bam, extra ->
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

    // Empirical strandedness detection
    ch_strandedness = INFER_STRANDEDNESS(ch_bams_strand, reference_gtf)
        .collect()

    // Per-sample StringTie2 assembly — collect GTFs for merge
    STRINGTIE2_PER_SAMPLE(ch_bams_assembly, reference_gtf)
    ch_gtf_list = STRINGTIE2_PER_SAMPLE.out.transcript_gtf
        .map { sample, gtf -> gtf }
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
        ch_cpat_models = Channel.value([
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
        STRINGTIE2_MERGE.out.merged_gtf
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
        ch_conservation = Channel.of("transcript_id\tmean_score\tmax_score\tpct_bases_conserved")
            .collectFile(name: 'conservation_placeholder.tsv', newLine: true).first()
    }

    if (params.synteny_chain_file && params.synteny_target_gtf) {
        SYNTENY_CHECK(ANNOTATION_FREEZE.out.frozen_gtf)
        ch_synteny = SYNTENY_CHECK.out.scores.first()
    } else {
        ch_synteny = Channel.of("transcript_id\tsyntenic_locus\tsyntenic_target_gene_id")
            .collectFile(name: 'synteny_placeholder.tsv', newLine: true).first()
    }

    // Gene catalog: one row per gene of the full analysis universe
    BUILD_GENE_CATALOG(
        ANNOTATION_FREEZE.out.frozen_gtf,
        reference_gtf,
        ch_conservation,
        ch_synteny
    )

    emit:
    frozen_gtf        = ANNOTATION_FREEZE.out.frozen_gtf
    tx2gene_detailed   = ANNOTATION_FREEZE.out.tx2gene_detailed
    gene_catalog       = BUILD_GENE_CATALOG.out.gene_catalog
    strand_consensus   = ch_strandedness
}
