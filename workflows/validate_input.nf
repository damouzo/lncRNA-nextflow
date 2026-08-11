// Input validation — samplesheet, contrast matrix, design formula sanity checks

workflow validateInput {
    take:
    input_csv
    genome_fa
    reference_gtf
    contrast_csv
    design_formula

    main:
    // Parse samplesheet: sample,condition,batch,fastq_1,fastq_2,bam,[covariates...]
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