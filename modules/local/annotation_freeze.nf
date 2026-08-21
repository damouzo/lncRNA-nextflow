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

    # Write frozen transcript ID list for GTF extraction
    nc_ids <- unique(nc\$ID)
    writeLines(nc_ids, "frozen_ids.txt")

    cat("Frozen lncRNAs:", length(nc_ids), "\\n")

    # Extract matching transcripts from source GTF by frozen transcript ID
    system(paste("grep -Ff frozen_ids.txt", shQuote("${source_gtf}"),
                 "> frozen_lncrna.gtf"))

    # Fail loudly instead of silently substituting the full candidate set
    if (file.info("frozen_lncrna.gtf")\$size < 100) {
        stop("No frozen transcripts matched in source GTF. Check transcript ID consistency between CPAT/CPC2 and the assembly.")
    }

    # Derive tx2gene directly from the frozen GTF so transcript/gene IDs match
    # exactly what gffread emits when building the Salmon index.
    gtf_lines <- readLines("frozen_lncrna.gtf")
    has_tx   <- grepl('transcript_id "', gtf_lines)
    tx_id    <- sub('.*transcript_id "([^"]+)".*', '\\\\1', gtf_lines[has_tx])
    gene_id  <- sub('.*gene_id "([^"]+)".*', '\\\\1', gtf_lines[has_tx])
    tx2gene  <- unique(data.frame(transcript_id = tx_id, gene_id = gene_id,
                                  stringsAsFactors = FALSE))
    tx2gene\$category <- "novel_lncRNA"

    write_tsv(tx2gene, "tx2gene_detailed.tsv")

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
