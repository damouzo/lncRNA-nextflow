// A10: Orthogonal catalog evidence overlap
// Tags transcripts: known_catalog_hit vs putatively_novel
// LNCipedia implemented; NONCODE/RNAcentral/FANTOM CAT stubbed for future extension
// Reporting only, never a filter
//
// Pure R implementation: rtracklayer + GenomicRanges handle GTF import and
// overlap detection. No external tools (gffread/bedtools) needed — the
// r_bioc container already ships both packages.

process CATALOG_OVERLAP {

    input:
    path consensus_table
    path consensus_gtf

    output:
    path "catalog_overlap.tsv", emit: overlap_table

    script:
    """
    #!/usr/bin/env Rscript
    library(rtracklayer)
    library(GenomicRanges)
    library(dplyr)
    library(readr)

    consensus <- read_tsv("${consensus_table}", show_col_types = FALSE)

    id_col <- names(consensus)[grepl("^ID\$|transcript_id", names(consensus), ignore.case = TRUE)][1]
    if (is.na(id_col)) id_col <- names(consensus)[1]

    consensus\$known_catalog_hit <- FALSE
    consensus\$catalog_sources   <- "none"

    # Import consensus GTF (transcript-only from CODING_CONSENSUS)
    gtf <- import("${consensus_gtf}")
    cat(sprintf("Imported %d transcript features from GTF\\n", length(gtf)))

    # Resolve transcript IDs from the imported GRanges: try mcols column
    # first, then row names as fallback
    tx_ids <- if ("transcript_id" %in% colnames(mcols(gtf))) {
        mcols(gtf)[["transcript_id"]]
    } else if (length(names(gtf)) > 0 && !all(grepl("^\\\\d+\$", names(gtf)))) {
        names(gtf)
    } else {
        character(0)
    }

    # ---- LNCipedia ----
    lnc_bed <- "${params.lncipedia_bed ?: ''}"
    if (nchar(lnc_bed) > 0 && lnc_bed != "null" && lnc_bed != "[]" && file.exists(lnc_bed)) {
        lncrna <- import(lnc_bed)
        ov <- findOverlaps(gtf, lncrna)
        ov_ids <- unique(tx_ids[queryHits(ov)])
        ov_ids <- ov_ids[ov_ids %in% consensus[[id_col]]]
        cat(sprintf("  LNCipedia: %d overlapping transcripts\\n", length(ov_ids)))
        if (length(ov_ids) > 0) {
            hit <- consensus[[id_col]] %in% ov_ids
            consensus\$known_catalog_hit[hit] <- TRUE
            consensus\$catalog_sources[hit]   <- "LNCipedia"
        }
    } else {
        cat("  LNCipedia BED not configured or not found — skipping\\n")
    }

    # ---- NONCODE (stub) ----
    # ---- RNAcentral (stub) ----
    # ---- FANTOM CAT (stub) ----

    cat(sprintf("known_catalog_hit TRUE: %d / %d\\n",
                sum(consensus\$known_catalog_hit), nrow(consensus)))
    write_tsv(consensus, "catalog_overlap.tsv")
    """
}