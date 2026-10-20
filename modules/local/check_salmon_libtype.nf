// B2b: Post-Salmon library-type guard (fail-fast, per sample)
//
// Hard criterion (when lib_format_counts.json lists per-format counts):
//   dominant format == declared libType
//   AND declared / (declared + opposite) >= salmon_min_orientation_ratio
//   (e.g. reverse paired-end: ISR / (ISR + ISF) >= 0.70)
// If per-format counts are not listed, that criterion degrades to a warning.
//
// compatible_fragment_ratio is informational (warning) with a hard floor
// salmon_min_compat_ratio: below the floor the run fails.
//
// errorStrategy 'terminate' aborts the cohort on the first bad sample.
// Regression case: the archived pre-fix quant dirs (declared SR, effective IU)
// must fail here.

process CHECK_SALMON_LIBTYPE {
    tag "${sample}"
    errorStrategy 'terminate'

    input:
    tuple val(sample), path(quant_dir)

    output:
    tuple val(sample), path("${sample}_quant"), emit: checked

    script:
    def min_ratio = (params.salmon_min_compat_ratio as double)
    def min_orient = (params.salmon_min_orientation_ratio as double)
    """
    #!/usr/bin/env bash
    set -euo pipefail
    python3 - "${sample}_quant" "${min_ratio}" "${min_orient}" <<'PYEOF'
    import json
    import os
    import sys

    qdir, min_ratio, min_orient = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
    sample = os.path.basename(qdir)
    LIB_KEYS = ["MSF", "OSF", "ISF", "MSR", "OSR", "ISR", "SF", "SR", "MU", "OU", "IU", "U"]

    def load(path):
        with open(path) as fh:
            return json.load(fh)

    cmd_path = os.path.join(qdir, "cmd_info.json")
    lfc_path = os.path.join(qdir, "lib_format_counts.json")
    meta_path = next((p for p in [os.path.join(qdir, "aux_info", "meta_info.json"),
                                  os.path.join(qdir, "meta_info.json")]
                      if os.path.exists(p)), None)
    for p in (cmd_path, lfc_path):
        if not os.path.exists(p):
            raise SystemExit("CHECK_SALMON_LIBTYPE [%s]: missing %s" % (sample, p))
    if meta_path is None:
        raise SystemExit("CHECK_SALMON_LIBTYPE [%s]: meta_info.json not found" % sample)

    cmd, lfc, meta = load(cmd_path), load(lfc_path), load(meta_path)
    declared = cmd.get("libType")
    expected = lfc.get("expected_format")
    ratio = lfc.get("compatible_fragment_ratio")
    counts = {k: int(lfc[k]) for k in LIB_KEYS if k in lfc}

    def warn(msg):
        print("WARNING CHECK_SALMON_LIBTYPE [%s]: %s" % (sample, msg))

    # ---- hard format checks (F0.9): a wrong -l aborts the cohort ----
    libtypes = meta.get("library_types") or []
    if expected != declared:
        raise SystemExit(
            "CHECK_SALMON_LIBTYPE [%s]: declared libType '%s' but expected_format is "
            "'%s'. Wrong library type." % (sample, declared, expected))
    if len(libtypes) != 1 or libtypes[0] != declared:
        raise SystemExit(
            "CHECK_SALMON_LIBTYPE [%s]: library_types %s inconsistent with declared "
            "'%s'." % (sample, libtypes, declared))

    # ---- dominant format (independent of the orientation ratio) ----
    if not counts:
        warn("per-format counts not listed in lib_format_counts.json; dominant-format criterion NOT enforced")
    elif declared in ("A", "", None) or not any(c in declared for c in "RFU"):
        warn("libType '%s' is auto/unknown; dominant-format criterion NOT enforced" % declared)
    else:
        opp = declared.replace("R", "F") if "R" in declared else (
            declared.replace("F", "R") if "F" in declared else None)
        dom = max(counts, key=counts.get)
        # If the declared format is not dominant the -l is wrong regardless of the
        # opposite-strand ratio (e.g. declared SR on paired data: SR count is 0).
        if dom != declared:
            raise SystemExit(
                "CHECK_SALMON_LIBTYPE [%s]: dominant format '%s' != declared '%s' "
                "(counts: %s). Wrong library type." % (sample, dom, declared, counts))
        if opp and opp in counts:
            d = counts.get(declared, 0)
            o = counts[opp]
            if (d + o) > 0:
                orient_ratio = d / (d + o)
                if orient_ratio < min_orient:
                    raise SystemExit(
                        "CHECK_SALMON_LIBTYPE [%s]: %s/(%s+%s)=%.3f < %.2f. Wrong or contaminated library."
                        % (sample, declared, declared, opp, orient_ratio, min_orient))
                print("CHECK_SALMON_LIBTYPE [%s]: dominant=%s %s/(%s+%s)=%.3f OK"
                      % (sample, dom, declared, declared, opp, orient_ratio))
        else:
            warn("counts for '%s' not listed; orientation-ratio criterion NOT enforced" % opp)

    # ---- compatible_fragment_ratio: warning with a hard floor ----
    if ratio is None:
        warn("compatible_fragment_ratio missing")
    elif float(ratio) < min_ratio:
        raise SystemExit(
            "CHECK_SALMON_LIBTYPE [%s]: compatible_fragment_ratio %.3f below hard floor %.2f"
            % (sample, float(ratio), min_ratio))
    elif float(ratio) < 1.0:
        warn("compatible_fragment_ratio %.3f < 1.0 (informational)" % float(ratio))
    PYEOF
    """
}
