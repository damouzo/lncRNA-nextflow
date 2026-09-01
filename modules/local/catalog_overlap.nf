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

    # ---- NONCODE ----
    ncode_bed <- "${params.noncode_bed ?: ''}"
    if (nchar(ncode_bed) > 0 && ncode_bed != "null" && ncode_bed != "[]" && file.exists(ncode_bed)) {
        ncode <- import(ncode_bed)
        ncode_ov <- findOverlaps(gtf, ncode)
        ncode_ov_ids <- unique(tx_ids[queryHits(ncode_ov)])
        ncode_ov_ids <- ncode_ov_ids[ncode_ov_ids %in% consensus[[id_col]]]
        cat(sprintf("  NONCODE: %d overlapping transcripts\\n", length(ncode_ov_ids)))
        if (length(ncode_ov_ids) > 0) {
            ncode_hit <- consensus[[id_col]] %in% ncode_ov_ids
            consensus\$known_catalog_hit[ncode_hit] <- TRUE
            consensus\$catalog_sources[ncode_hit] <- ifelse(
                consensus\$catalog_sources[ncode_hit] == "none",
                "NONCODE",
                paste(consensus\$catalog_sources[ncode_hit], "NONCODE", sep = ",")
            )
        }
    } else {
        cat("  NONCODE BED not configured or not found — skipping\\n")
    }

    # ---- RNAcentral (stub) ----
    # ---- FANTOM CAT (stub) ----

    cat(sprintf("known_catalog_hit TRUE: %d / %d\\n",
                sum(consensus\$known_catalog_hit), nrow(consensus)))
    write_tsv(consensus, "catalog_overlap.tsv")
    """
}