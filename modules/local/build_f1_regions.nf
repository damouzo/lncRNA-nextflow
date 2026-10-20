// F1.0: build candidate + control regions for the F1.3/F1.4 metrics.
// Phase A, post-freeze, terminal (nothing depends on its output). Pure Python
// interval logic (bin/f1_regions.py). Deterministic (fixed seed).

process BUILD_F1_REGIONS {

    input:
    path frozen_gtf
    path reference_gtf
    path recurrence_table

    output:
    path "f1_regions/regions.tsv", emit: regions
    path "f1_regions/flank.saf", emit: flank_saf
    path "f1_regions/antisense.saf", emit: antisense_saf
    path "f1_regions/regions_stats.tsv", emit: regions_stats

    script:
    """
    f1_regions.py \\
        --frozen-gtf "${frozen_gtf}" \\
        --reference-gtf "${reference_gtf}" \\
        --recurrence "${recurrence_table}" \\
        --outdir f1_regions \\
        --flank-bp ${params.f1_flank_bp} \\
        --min-flank-len ${params.f1_min_flank_len} \\
        --cap ${params.f1_control_cap} \\
        --seed ${params.f1_seed} \\
        --expr-min-reps ${params.expression_recurrence_min_reps} \\
        --tier1-classes '${params.tier1_classes.join(",")}' \\
        --tier1a-classes '${params.tier1a_classes.join(",")}' \\
        --tier1-parent-biotypes '${params.tier1_parent_biotypes.join(",")}' \\
        --annotated-lncrna-biotypes '${params.annotated_lncrna_biotypes.join(",")}'
    """
}
