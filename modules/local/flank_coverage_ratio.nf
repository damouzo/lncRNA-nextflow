// F1.3b: aggregate per-sample flank counts and compute the flank metrics.
// One task for the whole cohort: merges the collected single-column matrices and
// derives per-element, per-norm_group metrics. Metric only; thresholds in
// DERIVE_THRESHOLDS so a rule change re-runs from -resume without re-counting.

process FLANK_COVERAGE_RATIO {

    input:
    path s2_files             // collected per-sample -s 2 tables
    path s1_files             // collected per-sample -s 1 tables
    path coldata              // sample, norm_group, bam basename
    path regions

    output:
    path "flank_metrics.tsv", emit: flank_metrics

    script:
    """
    f1_metrics.py flank \\
        --s2 ${s2_files.join(' ')} \\
        --s1 ${s1_files.join(' ')} \\
        --regions "${regions}" \\
        --coldata "${coldata}" \\
        --out flank_metrics.tsv \\
        --min-reads ${params.f1_min_reads}
    """
}
