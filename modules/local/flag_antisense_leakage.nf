// F1.4b: aggregate per-sample antisense counts and compute the leakage metrics.
// One task for the whole cohort; metric only (thresholds in DERIVE_THRESHOLDS).

process FLAG_ANTISENSE_LEAKAGE {

    input:
    path s2_files             // collected per-sample -s 2 tables
    path s1_files             // collected per-sample -s 1 tables
    path coldata              // sample, norm_group, bam basename
    path regions

    output:
    path "antisense_metrics.tsv", emit: antisense_metrics

    script:
    """
    f1_metrics.py antisense \\
        --s2 ${s2_files.join(' ')} \\
        --s1 ${s1_files.join(' ')} \\
        --regions "${regions}" \\
        --coldata "${coldata}" \\
        --out antisense_metrics.tsv \\
        --min-reads ${params.f1_min_reads}
    """
}
