// F6.1 / F6.2: reviewer validations from the catalog + pooled expression.
// Cheap, runs once at the end of the pipeline; publishes to deliverables/validation/.

process RECOVERY_CHECK {

    label 'process_single'

    input:
    path catalog_genes      // deliverables catalog_genes.tsv
    path tpm_tsvs           // collected deliverables/expression/<group>_tpm.tsv

    output:
    path "recovery_reference.tsv", emit: recovery
    path "tier1_vs_annotated.tsv", emit: comparison
    path "versions.yml", emit: versions

    script:
    def tpm_args = tpm_tsvs.collect { "--tpm \"${it}\"" }.join(' ')
    """
    recovery_analysis.py --catalog "${catalog_genes}" ${tpm_args} --outdir .

    emit_versions.sh "${task.process}" "${task.container}" python3 python3 --version
    """
}
