// Input validation — samplesheet, contrast matrix, design formula sanity checks
// Supports multi-lane FASTQs: one row per sample-lane, grouped by sample downstream

workflow validateInput {
    take:
    input_csv
    genome_fa
    reference_gtf
    contrast_csv
    design_formula

    main:
    // Parse samplesheet: sample,condition,batch,fastq_1,fastq_2,bam,[covariates...]
    // One row per SAMPLE-LANE — multi-lane samples have multiple rows with same sample ID
    ch_samplesheet = Channel.fromPath(input_csv, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def extra = [:]
            row.each { col, val ->
                if (!['sample', 'condition', 'batch', 'fastq_1', 'fastq_2', 'bam'].contains(col)) {
                    extra[col] = val
                }
            }
            tuple(
                row.sample,
                row.condition,
                row.batch,
                file(row.fastq_1, checkIfExists: true),
                file(row.fastq_2, checkIfExists: true),
                file(row.bam, checkIfExists: true),
                extra
            )
        }
        // Group lanes by sample: FASTQ lists aggregated, metadata taken from first lane
        .groupTuple(by: 0)
        .map { sample, conds, batches, fq1s, fq2s, bams, extras ->
            tuple(
                sample,
                conds[0],
                batches[0],
                fq1s,       // List[Path] — multiple lanes
                fq2s,       // List[Path]
                bams[0],    // BAM is same for all lanes
                extras[0]
            )
        }
        .set { ch_parsed }

    // Parse contrast matrix: contrast_name, numerator, denominator, batch(optional)
    ch_contrasts = Channel.fromPath(contrast_csv, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def batch = row.containsKey('batch') ? (row.batch ?: null) : null
            tuple(row.contrast_name, row.numerator, row.denominator, batch)
        }
        .set { ch_contrasts }

    // Validate required params
    if (!design_formula || design_formula.trim().isEmpty()) {
        error("design_formula is required (e.g. '~ batch + condition')")
    }
    if (!genome_fa) { error("genome FASTA is required") }
    if (!reference_gtf) { error("reference GTF is required") }

    emit:
    ch_samplesheet  = ch_parsed
    ch_contrasts    = ch_contrasts
    design_formula   = design_formula
    coldata_csv      = file(input_csv, checkIfExists: true)
    contrasts_csv    = file(contrast_csv, checkIfExists: true)
}