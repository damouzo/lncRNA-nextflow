// Split samplesheet and contrasts by norm_group
// Writes per-group coldata CSV files (deduplicated, one row per sample) for Phase B

process SPLIT_BY_NORMGROUP {
    label 'process_single'

    input:
    path samplesheet_csv
    path contrasts_csv

    output:
    path "coldata_*.csv", emit: per_group

    script:
    """
    #!/usr/bin/env Rscript
    library(readr)

    coldata <- read.csv("${samplesheet_csv}", stringsAsFactors = FALSE, check.names = FALSE)
    stopifnot("norm_group" %in% colnames(coldata),
        "samplesheet must have a 'norm_group' column")

    groups <- sort(unique(coldata[["norm_group"]]))
    for (g in groups) {
        sub <- coldata[coldata[["norm_group"]] == g, , drop = FALSE]
        # Drop lane-level columns and deduplicate by sample
        drop_cols <- intersect(c("fastq_1", "fastq_2", "bam"), colnames(sub))
        sub <- sub[!duplicated(sub[["sample"]]),
                   setdiff(colnames(sub), drop_cols),
                   drop = FALSE]
        outfile <- paste0("coldata_", g, ".csv")
        write.csv(sub, outfile, row.names = FALSE, quote = FALSE)
        cat(sprintf("Wrote %s (%d samples)\\n", outfile, nrow(sub)))
    }
    """
}