// Slice D: optional xlsx export of the per-contrast DE table.
// errorStrategy 'ignore': a formatting/Excel failure must never fail the run.
// Gated by params.export_xlsx upstream. The TSV remains the source of truth.

process MAKE_XLSX {

    label 'process_single'
    errorStrategy 'ignore'

    input:
    tuple val(group), val(ct), path(de_tsv)

    output:
    path "${group}__${ct}.xlsx", emit: xlsx

    script:
    """
    make_xlsx.R --in "${de_tsv}" --out "${group}__${ct}.xlsx" --sheet "${ct}"
    """
}
