// A1: BAM content validation via samtools quickcheck
// Hard-fail if any BAM is corrupt or truncated

process BAM_QUICKCHECK {
    tag "${sample}"

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path(bam), emit: validated
    path "${sample}.quickcheck.log", emit: logs

    script:
    """
    samtools quickcheck -v "${bam}" 2>&1 | tee "${sample}.quickcheck.log"
    if grep -q "ERROR" "${sample}.quickcheck.log"; then
        echo "FATAL: BAM validation failed for ${sample}" >&2
        exit 1
    fi
    """
}
