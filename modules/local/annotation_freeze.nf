// Annotation freezing — frozen novel GTF + full-universe tx2gene + checksum manifest
// Strict barrier: Phase B only receives frozen artifacts.
// tx2gene_detailed.tsv now spans the entire quantification universe
// (all reference transcripts + novel lncRNAs) so DESeq2 normalizes and
// corrects over the full transcriptome, matching the Salmon index.

process ANNOTATION_FREEZE {

    input:
    path recurrence_table
    path source_gtf            // length-filtered GTF with all candidate transcripts
    path genome_fa
    path reference_gtf         // reference annotation (full tx2gene universe)

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
    system(paste("grep -F -f frozen_ids.txt", shQuote("${source_gtf}"),
                 "> frozen_lncrna.gtf"))

    # Fail loudly instead of silently substituting the full candidate set
    if (file.info("frozen_lncrna.gtf")\$size < 100) {
        stop("No frozen transcripts matched in source GTF. Check transcript ID consistency between CPAT/CPC2 and the assembly.")
    }

    # Keep these two helpers in sync with BUILD_GENE_CATALOG: both must strip
    # Ensembl versions but leave MSTRG IDs verbatim, or downstream gene_id
    # joins (tx2gene / catalog / DE / report) break silently.
    extract_attr <- function(attrs, key) {
        pat <- paste0('(?:^|; )', key, ' "([^"]+)"')
        m <- regexpr(pat, attrs, perl = TRUE)
        ans <- rep(NA_character_, length(attrs))
        ok <- m > 0L
        if (any(ok)) ans[ok] <- sub(paste0('.*', pat, '.*\$'), '\\\\1', attrs[ok], perl = TRUE)
        ans
    }
    strip_version <- function(x) sub("\\\\.[0-9]+\$", "", x)

    # ---- Reference tx2gene: one transcript per line, same universe as the index ----
    system(paste("grep -P '\\ttranscript\\t'", shQuote("${reference_gtf}"),
                 "> ref_transcript_lines.txt"))
    ref_lines <- readLines("ref_transcript_lines.txt")
    ref_parts <- strsplit(ref_lines, "\\t", fixed = TRUE)
    ref_ok    <- lengths(ref_parts) >= 9L
    ref_attrs <- vapply(ref_parts[ref_ok], function(p) p[[9L]], character(1))
    ref_txg <- unique(data.frame(
        transcript_id = strip_version(extract_attr(ref_attrs, "transcript_id")),
        gene_id       = strip_version(extract_attr(ref_attrs, "gene_id")),
        stringsAsFactors = FALSE
    ))
    ref_txg <- ref_txg[!is.na(ref_txg\$transcript_id) & !is.na(ref_txg\$gene_id), ]
    ref_txg\$category <- "reference"

    # ---- Novel tx2gene from the frozen GTF ----
    # MSTRG IDs are not versioned like Ensembl — keep them verbatim
    gtf_lines <- readLines("frozen_lncrna.gtf")
    gtf_parts <- strsplit(gtf_lines, "\\t", fixed = TRUE)
    gtf_ok    <- lengths(gtf_parts) >= 9L

    # `stringtie --merge -G` reuses the reference transcript_id for transcripts
    # that only partially match it (class_code != "=", e.g. containment "c"/"k"),
    # which then collides with the same ID in ref_txg above. tximport keeps only
    # the first tx2gene match per transcript, so every read for that transcript
    # silently gets attributed to the reference gene and the novel gene ends up
    # with zero counts everywhere. Disambiguate those IDs before they're used
    # for quantification or tx2gene; exact matches ("=") keep the real accession.
    ok_lines  <- gtf_lines[gtf_ok]
    ok_parts  <- strsplit(ok_lines, "\\t", fixed = TRUE)
    feature   <- vapply(ok_parts, function(p) p[[3L]], character(1))
    attrs_all <- vapply(ok_parts, function(p) p[[9L]], character(1))

    tx_mask    <- feature == "transcript"
    tx_ids     <- extract_attr(attrs_all[tx_mask], "transcript_id")
    tx_genes   <- extract_attr(attrs_all[tx_mask], "gene_id")
    tx_class   <- extract_attr(attrs_all[tx_mask], "class_code")
    ref_id_set <- unique(ref_txg\$transcript_id)

    needs_rename <- !is.na(tx_class) & tx_class != "=" & strip_version(tx_ids) %in% ref_id_set
    if (any(needs_rename)) {
        old_ids    <- tx_ids[needs_rename]
        rename_map <- setNames(paste0(tx_genes[needs_rename], "__", old_ids), old_ids)

        line_ids <- extract_attr(attrs_all, "transcript_id")
        hit      <- line_ids %in% old_ids
        ok_lines[hit] <- mapply(function(line, old) {
            sub(paste0('transcript_id "', old, '"'),
                paste0('transcript_id "', rename_map[[old]], '"'),
                line, fixed = TRUE)
        }, ok_lines[hit], line_ids[hit])

        cat(sprintf("INFO: disambiguated %d transcript IDs reused from the reference (class_code != '=')\\n",
                    length(old_ids)))
        gtf_lines[gtf_ok] <- ok_lines
        writeLines(gtf_lines, "frozen_lncrna.gtf")
    }

    gtf_attrs <- vapply(strsplit(ok_lines, "\\t", fixed = TRUE), function(p) p[[9L]], character(1))
    has_tx <- grepl("transcript_id", gtf_attrs, fixed = TRUE)
    gtf_attrs <- gtf_attrs[has_tx]
    nv_tx2gene <- unique(data.frame(
        transcript_id = extract_attr(gtf_attrs, "transcript_id"),
        gene_id       = extract_attr(gtf_attrs, "gene_id"),
        stringsAsFactors = FALSE
    ))
    nv_tx2gene <- nv_tx2gene[!is.na(nv_tx2gene\$transcript_id) & !is.na(nv_tx2gene\$gene_id), ]
    nv_tx2gene\$category <- "novel_lncRNA"

    tx2gene <- unique(rbind(ref_txg, nv_tx2gene))

    write_tsv(tx2gene, "tx2gene_detailed.tsv")

    # Checksums
    sha_gtf <- system("sha256sum frozen_lncrna.gtf | cut -d' ' -f1", intern = TRUE)
    sha_tsv <- system("sha256sum tx2gene_detailed.tsv | cut -d' ' -f1", intern = TRUE)

    sink("frozen_manifest.json")
    cat(sprintf('{
  "pipeline_version": "0.1.0",
  "freeze_date": "%s",
  "n_novel_lncrna": %d,
  "n_reference_transcripts": %d,
  "n_total_transcripts": %d,
  "sha256": {
    "frozen_gtf": "%s",
    "tx2gene_detailed": "%s"
  }
}
', Sys.time(), nrow(nc), nrow(ref_txg), nrow(tx2gene),
         trimws(sha_gtf), trimws(sha_tsv)))
    sink()

    writeLines(c(trimws(sha_gtf), trimws(sha_tsv)), "frozen_checksums.sha256")
    """
}