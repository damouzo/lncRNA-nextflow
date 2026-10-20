// A2b: Strandedness input guard (fail-fast)
//
// Compares the per-sample RSeQC outputs against the declared library_strandedness
// contract BEFORE Salmon runs, and reports the normalised fidelity E/(E+O).
//
// Thin wrapper: all logic lives in bin/check_strandedness.py (unit-tested with
// synthetic RSeQC outputs in bin/tests/test_strandedness.py). This is the
// ORIENTATION guard (check_salmon_libtype only echoes the -l it was given).
//
// Fails when:
//   - samples disagree on paired/single-end library type
//   - the majority empirical orientation does not match library_strandedness
//     (declaring forward / unstranded on reverse data stops the run)
//   - strandedness_guard_strict=true and any sample is below the fidelity floor
// Warns (does not fail) when a sample is below strandedness_min_fidelity.
//
// Writes run_info/strandedness_summary.tsv. Never touches frozen_manifest.json.

process CHECK_STRANDEDNESS {

    input:
    path strandedness_files

    output:
    path "strandedness_summary.tsv", emit: summary

    script:
    def declared = Strandedness.mode(params.library_strandedness)
    def min_fid  = (params.strandedness_min_fidelity as double)
    def strict   = params.strandedness_guard_strict ? '--strict' : ''
    """
    check_strandedness.py --files *.strandedness.txt --declared ${declared} \
        --min-fidelity ${min_fid} ${strict} --summary-out strandedness_summary.tsv
    """
}
