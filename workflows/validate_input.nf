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
            // Every lane of a sample must agree on condition and batch: a mismatch
            // means the samplesheet is wrong and the design would silently use lane 1.
            if (conds.unique().size() > 1) {
                error("sample '${sample}': condition differs across lanes (${conds.unique().join(', ')})")
            }
            if (batches.unique().size() > 1) {
                error("sample '${sample}': batch differs across lanes (${batches.unique().join(', ')})")
            }
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
    // The numerator/denominator are checked against the conditions actually present
    // in the samplesheet, so a typo fails before any quantification runs.
    // wrap the collected list so combine() keeps it as one field; otherwise it
    // flattens the conditions into the contrast tuple and the map below breaks
    ch_conditions = ch_parsed
        .map { _sample, condition, _batch, _fq1, _fq2, _bam, _extra -> condition }
        .unique()
        .collect()
        .map { conds -> [ conds ] }
    ch_contrasts = Channel.fromPath(contrast_csv, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def batch = row.containsKey('batch') ? (row.batch ?: null) : null
            tuple(row.contrast_name, row.numerator, row.denominator, batch)
        }
        .combine(ch_conditions)
        .map { name, numerator, denominator, batch, conditions ->
            if (!(numerator in conditions)) {
                error("contrast '${name}': numerator '${numerator}' is not a condition in the samplesheet " +
                      "(known: ${conditions.join(', ')})")
            }
            if (!(denominator in conditions)) {
                error("contrast '${name}': denominator '${denominator}' is not a condition in the samplesheet " +
                      "(known: ${conditions.join(', ')})")
            }
            tuple(name, numerator, denominator, batch)
        }

    // Validate required params
    if (!design_formula || design_formula.trim().isEmpty()) {
        error("design_formula is required (e.g. '~ batch + condition')")
    }
    if (!genome_fa) { error("genome FASTA is required") }
    if (!reference_gtf) { error("reference GTF is required") }

    // The freeze gate (novel_classes) and the Tier definition are separate params,
    // but a Tier 1 class that is never frozen could not exist in the catalog.
    def novel_cls = params.novel_classes as List
    def tier1_outside = (params.tier1_classes as List).findAll { !(it in novel_cls) }
    if (!novel_cls) { error("novel_classes must not be empty") }
    if (tier1_outside) {
        error("tier1_classes ${tier1_outside} are not in novel_classes ${novel_cls}: they would never be frozen")
    }

    emit:
    ch_samplesheet  = ch_parsed
    ch_contrasts    = ch_contrasts
    design_formula   = design_formula
    coldata_csv      = file(input_csv, checkIfExists: true)
    contrasts_csv    = file(contrast_csv, checkIfExists: true)
}