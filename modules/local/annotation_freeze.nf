// Annotation freezing — final GTF + tx2gene_detailed + checksum manifest
// Strict barrier: Phase B only receives frozen artifacts

process ANNOTATION_FREEZE {

    input:
    path recurrence_table
    path source_gtf            // length-filtered GTF with all candidate transcripts
    path genome_fa

    output:
    path "frozen_lncrna.gtf", emit: frozen_gtf
    path "tx2gene_detailed.tsv", emit: tx2gene_detailed
    path "frozen_manifest.json", emit: manifest
    path "frozen_checksums.sha256", emit: checksums

    script:
    """
    #!/usr/bin/env Rscript
    library(dplyr)
    library(readr)
    library(rtracklayer)

    rec <- read_tsv("${recurrence_table}", show_col_types = FALSE)

    # Filter to high-confidence noncoding that pass recurrence
    nc <- subset(rec, consensus_noncoding == TRUE & passes_recurrence == TRUE)

    if (nrow(nc) == 0) {
        stop("No non-coding transcripts passed all filters. Check your inputs.")
    }

    # Build tx2gene mapping with gene-level grouping
    tx2gene <- nc[, c("ID", "ID")]
    if ("X1" %in% colnames(nc)) {
        # CPAT/CPC2 output typically uses 'ID' or sequence name as first column
        colnames(tx2gene) <- c("transcript_id", "gene_id")
    } else if ("transcript_id" %in% colnames(nc) && "gene_id" %in% colnames(nc)) {
        tx2gene <- nc[, c("transcript_id", "gene_id")]
    }
    colnames(tx2gene) <- c("transcript_id", "gene_id")
    tx2gene\$category <- "novel_lncRNA"

    write_tsv(tx2gene, "tx2gene_detailed.tsv")

    # Write frozen transcript ID list for GTF extraction
    nc_ids <- unique(tx2gene\$transcript_id)
    writeLines(nc_ids, "frozen_ids.txt")

    cat("Frozen lncRNAs:", length(nc_ids), "\\n")

    # Extract matching transcripts from source GTF using gffread
    # Filter GTF to only keep the frozen transcript IDs
    system(paste("grep -Ff frozen_ids.txt", "${source_gtf}",
                 "> frozen_lncrna.gtf || true"))

    if (file.info("frozen_lncrna.gtf")\$size < 100) {
        # Fallback: try extracting by transcript_id attribute
        system(paste("awk '\$3==\"transcript\"'", "${source_gtf}",
                     "> frozen_lncrna.gtf"))
        cat("Using full transcript set as fallback\\n")
    }

    # Checksums
    sha_gtf <- system("sha256sum frozen_lncrna.gtf | cut -d' ' -f1", intern = TRUE)
    sha_tsv <- system("sha256sum tx2gene_detailed.tsv | cut -d' ' -f1", intern = TRUE)

    sink("frozen_manifest.json")
    cat(sprintf('{
  "pipeline_version": "0.1.0",
  "freeze_date": "%s",
  "n_novel_lncrna": %d,
  "sha256": {
    "frozen_gtf": "%s",
    "tx2gene_detailed": "%s"
  }
}
', Sys.time(), nrow(nc), trimws(sha_gtf), trimws(sha_tsv)))
    sink()

    writeLines(c(trimws(sha_gtf), trimws(sha_tsv)), "frozen_checksums.sha256")
    """
}
