// F1.5a: derive flag thresholds from the controls (cheap; reads only metric TSVs).
// Overrides win; unreliable controls -> config defaults + thresholds_unreliable.
// Terminal.

process DERIVE_THRESHOLDS {

    input:
    path flank_metrics
    path antisense_metrics

    output:
    path "thresholds_used.tsv", emit: thresholds
    path "sensitivity_thresholds.tsv", emit: sensitivity
    path "thresholds_used.txt", emit: thresholds_txt
    path "f1_controls", emit: controls

    script:
    def flank_override = params.flank_ratio_threshold != null ? "--flank-override ${params.flank_ratio_threshold}" : ""
    def anti_override  = params.antisense_ratio_threshold != null ? "--antisense-override ${params.antisense_ratio_threshold}" : ""
    def no_flank       = params.run_flank_check ? "" : "--no-flank"
    """
    f1_thresholds.py \\
        --flank-metrics "${flank_metrics}" \\
        --antisense-metrics "${antisense_metrics}" \\
        --outdir . \\
        --flank-default ${params.flank_ratio_default} \\
        --antisense-default ${params.antisense_ratio_default} \\
        --flank-percentile ${params.f1_flank_percentile} \\
        --antisense-percentile ${params.f1_antisense_percentile} \\
        --min-controls ${params.f1_min_controls} \\
        --neg-reject-min ${params.f1_neg_reject_min} \\
        --pos-flag-max ${params.f1_antisense_pos_flag_max} \\
        ${no_flank} ${flank_override} ${anti_override}
    """
}
