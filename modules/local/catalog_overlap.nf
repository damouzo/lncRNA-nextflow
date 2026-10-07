// A10: Orthogonal catalog evidence overlap (two levels)
// Tags transcripts with:
//   - STRICT level (only source of known_catalog_hit / novel figures): transcript
//     overlap on the CANDIDATE's strand, reciprocal threshold
//     external_catalog_min_reciprocal. Never filters.
//   - INFORMATIVE level (never filters): per-source bp coverage of the candidate
//     locus as a fraction of the candidate span.
// Strand is anchored to the CANDIDATE, never to the reference (for class codes
// x/n/k/s/e "same strand as candidate" == opposite of reference). bedtools
// -s is therefore NOT applied here.
// Sources: LNCipedia, NONCODE, RNAcentral, FANTOM CAT (any missing BED -> NA/0,
// never an error). RNAcentral is best-effort (expect a low hit rate).
//
// Pure R: rtracklayer + GenomicRanges (r_bioc container). No external binaries.

process CATALOG_OVERLAP {

    input:
    path consensus_table
    path consensus_gtf

    output:
    path "catalog_overlap.tsv", emit: overlap_table

    script:
    """
    #!/usr/bin/env Rscript
    suppressWarnings({
        library(rtracklayer)
        library(GenomicRanges)
        library(dplyr)
        library(readr)
    })

    consensus <- read_tsv("${consensus_table}", show_col_types = FALSE)

    id_col <- names(consensus)[grepl("^ID\$|transcript_id", names(consensus), ignore.case = TRUE)][1]
    if (is.na(id_col)) id_col <- names(consensus)[1]
    if (is.na(id_col) || !nzchar(id_col)) stop("CATALOG_OVERLAP: no transcript id column in the consensus table")

    consensus\$known_catalog_hit <- FALSE
    consensus\$catalog_sources   <- "none"
    for (sc in c("lncipedia", "noncode", "rnacentral", "fantomcat")) {
        consensus[[paste0(sc, "_id")]]  <- NA_character_
        consensus[[paste0(sc, "_pct")]] <- 0
    }

    # Import consensus GTF: standard GTF carries transcript AND exon features, so
    # keep transcript rows only — one feature per transcript for the overlap.
    gtf <- import("${consensus_gtf}")
    gtf <- gtf[mcols(gtf)\$type == "transcript"]
    gtf <- gtf[seqnames(gtf) %in% seqlevelsInUse(gtf) & width(ranges(gtf)) > 0]
    cat(sprintf("Imported %d transcript features from GTF\\n", length(gtf)))

    # Transcript id per GTF position. Hits are mapped back to the consensus
    # table BY ID, never by row position: the GTF keeps file order while the
    # consensus table is ID-sorted and may contain non-consensus rows too.
    tx_ids_pos <- if ("transcript_id" %in% colnames(mcols(gtf))) {
        as.character(mcols(gtf)[["transcript_id"]])
    } else if (length(names(gtf)) > 0 && !all(grepl("^\\\\d+\$", names(gtf)))) {
        names(gtf)
    } else {
        stop("CATALOG_OVERLAP: no transcript_id in the consensus GTF — cannot map hits")
    }
    if (length(tx_ids_pos) != length(gtf)) stop("CATALOG_OVERLAP: transcript_id length mismatch in consensus GTF")
    if (any(duplicated(tx_ids_pos))) stop("CATALOG_OVERLAP: duplicate transcript_id in consensus GTF")
    miss_ids <- tx_ids_pos[!tx_ids_pos %in% consensus[[id_col]]]
    if (length(miss_ids) > 0) {
        cat(sprintf("WARNING: %d consensus GTF ids missing from the consensus table\\n", length(miss_ids)))
    }

    MIN_STRICT <- as.numeric("${params.external_catalog_min_reciprocal}")

    # Annotate the consensus table from one external source BED.
    # 'strict' = same strand where both are known AND reciprocal coverage
    # >= MIN_STRICT (the parameter documented as the strict threshold);
    # raw coverage itself stays strand-agnostic.
    annotate_source <- function(ext, src) {
        ids <- rep(NA_character_, nrow(consensus))
        cov <- rep(0, nrow(consensus))
        coord_to <- match(tx_ids_pos, consensus[[id_col]])
        hits <- findOverlaps(gtf, ext, type = "any")
        if (length(hits) > 0) {
            q <- queryHits(hits); s <- subjectHits(hits)
            qs <- as.character(strand(gtf)[q]); ss <- as.character(strand(ext)[s])
            ext_id <- if ("name" %in% colnames(mcols(ext))) {
                z <- as.character(mcols(ext)\$name[s]); z[is.na(z)] <- src; z
            } else as.character(s)
            # ranges() drops strand metadata, so no ignore.strand arg (unknown in this
            # IRanges API version).
            int_w <- width(pintersect(ranges(gtf)[q], ranges(ext)[s]))
            span_w <- width(ranges(gtf))[q]
            for (k in seq_along(q)) {
                crow <- coord_to[q[k]]
                if (is.na(crow)) next
                cov[crow] <- min(1, cov[crow] + int_w[k] / span_w[k])
            }
            hit_row <- coord_to[q]
            sel <- which(!is.na(hit_row) &
                         (is.na(qs) | qs == "*" | is.na(ss) | ss == "*" | qs == ss) &
                         cov[hit_row] >= MIN_STRICT)
            if (length(sel) > 0) {
                bypos <- split(ext_id[sel], q[sel])
                for (qi in names(bypos)) {
                    crow <- coord_to[as.integer(qi)]
                    ids[crow] <- paste(unique(bypos[[qi]]), collapse = ",")
                }
            }
        }
        list(ids = ids, cov = cov)
    }

    sources <- list(
        lncipedia  = list(param = "${params.lncipedia_bed ?: ''}"),
        noncode    = list(param = "${params.noncode_bed ?: ''}"),
        rnacentral = list(param = "${params.rnacentral_bed ?: ''}"),
        fantomcat  = list(param = "${params.fantomcat_bed ?: ''}")
    )

    for (sc in names(sources)) {
        bed <- sources[[sc]]\$param
        if (nchar(bed) == 0 || bed %in% c("null", "[]", "none") || !file.exists(bed)) {
            cat(sprintf("  %s: BED not configured or not found — skipped\\n", sc))
            next
        }
        ext <- tryCatch(import(bed), error = function(e) NULL)
        if (is.null(ext)) {
            cat(sprintf("  %s: failed to import BED — skipped\\n", sc))
            next
        }
        res <- annotate_source(ext, sc)
        consensus[[paste0(sc, "_id")]]  <- res\$ids
        consensus[[paste0(sc, "_pct")]] <- res\$cov
        hit <- !is.na(res\$ids)
        consensus\$known_catalog_hit[hit] <- TRUE
        consensus\$catalog_sources[hit] <- ifelse(
            consensus\$catalog_sources[hit] == "none", sc,
            paste(consensus\$catalog_sources[hit], sc, sep = ","))
        cat(sprintf("  %s: %d strict transcripts, max coverage %.2f\\n",
                    sc, sum(hit, na.rm = TRUE), max(res\$cov, 0)))
    }

    cat(sprintf("known_catalog_hit TRUE: %d / %d\\n",
                sum(consensus\$known_catalog_hit), nrow(consensus)))
    write_tsv(consensus, "catalog_overlap.tsv")
    """
}