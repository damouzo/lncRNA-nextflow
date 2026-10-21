// Tool versions for the heavy Phase A/Fase B images, collected once per image.
// Deliberately kept OUT of the tool processes: editing SALMON_QUANTIFY,
// STRINGTIE2_*, GFFCOMPARE, FEATURECOUNTS_EXPRESSION, CPAT or CPC2 would
// invalidate their task cache and force a near-full recompute on -resume.
// One task per image, no data dependency, run once.

process VERSION_DISCOVERY {

    label 'process_single'

    output:
    path "versions.yml", emit: versions

    script:
    """
    emit_versions.sh "SALMON_QUANTIFY" "${task.container}" salmon salmon --version
    emit_versions.sh "STRINGTIE2_PER_SAMPLE" "${task.container}" stringtie stringtie --version
    emit_versions.sh "STRINGTIE2_MERGE" "${task.container}" stringtie stringtie --version
    emit_versions.sh "GFFCOMPARE" "${task.container}" gffcompare gffcompare --version
    emit_versions.sh "FEATURECOUNTS_EXPRESSION" "${task.container}" featureCounts featureCounts -v
    """
}

process VERSION_CODING {

    label 'process_single'

    output:
    path "versions.yml", emit: versions

    script:
    """
    emit_versions.sh "CPAT" "${task.container}" cpat.py cpat.py --version
    emit_versions.sh "CPC2" "${task.container}" CPC2.py CPC2.py --version
    """
}
