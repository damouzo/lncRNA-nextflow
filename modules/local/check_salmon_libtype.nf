// B2b: Post-Salmon library-type guard (fail-fast, per sample)
//
// Thin wrapper: all logic lives in bin/check_salmon_libtype.py (unit-tested with
// synthetic Salmon JSON in bin/tests/test_salmon_guard.py). Hard criteria:
//   H1  expected_format == cmd_info.libType   (detects Salmon ignoring the -l)
//   H2  meta_info.library_types == [libType]
//   H3  end-type coherence (2 read files <-> paired libtype)
// Reported as QC only: consistent_fraction, compatible_fragment_ratio, dominant
// format. Orientation is decided by CHECK_STRANDEDNESS (RSeQC), not here.
//
// errorStrategy 'terminate' aborts the cohort on the first bad sample (the
// regression case: the archived pre-fix quant dirs, declared SR on PE data).

process CHECK_SALMON_LIBTYPE {
    tag "${sample}"
    errorStrategy 'terminate'

    input:
    tuple val(sample), path(quant_dir)

    output:
    tuple val(sample), path("${sample}_quant"), emit: checked

    script:
    """
    check_salmon_libtype.py --quant ${sample}_quant
    """
}
