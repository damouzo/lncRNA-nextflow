// Annotation freezing — frozen novel catalog GTF + quantification universe + manifest
// Strict barrier: Phase B only receives frozen artifacts.
//
// Two GTFs come out of here, and they differ in exactly one way:
//   frozen_lncrna.gtf    the CATALOG entity set: novel_classes minus isoform_classes.
//                        Read by every catalog consumer (registry, gene catalog,
//                        F1/tier, publish). A reference-derived novel isoform (`j`)
//                        is not an independent locus, so it is absent here.
//   quant_universe.gtf   the QUANTIFICATION set: the catalog set plus the isoforms.
//                        Read only by BUILD_DECOY_INDEX, so transcript-level EM still
//                        distributes reads between the reference isoforms and the
//                        novel isoform. BUILD_TX2GENE remaps those isoform transcripts
//                        to their reference parent gene, so the parent total is not
//                        diluted and no isoform is ever counted as its own "lncRNA".
//
// tx2gene is NOT built here: BUILD_TX2GENE is its single producer.

process ANNOTATION_FREEZE {

    input:
    path recurrence_table
    path source_gtf            // length-filtered GTF with all candidate transcripts
    path genome_fa
    path reference_gtf         // reference annotation (kept for the manifest provenance)

    output:
    path "frozen_lncrna.gtf",  emit: frozen_gtf         // catalog set (no isoforms)
    path "quant_universe.gtf", emit: quant_universe     // catalog + isoforms (index set)
    path "frozen_manifest.json", emit: manifest
    path "frozen_checksums.sha256", emit: checksums

    script:
    """
    #!/usr/bin/env Rscript
    library(dplyr)
    library(readr)
    library(rtracklayer)

    rec <- read_tsv("${recurrence_table}", show_col_types = FALSE)

    # Freeze only the classes listed in novel_classes. The reference-derived
    # classes (c/k/m/n, and =, e, p, ...) are variants of an annotated transcript:
    # freezing them would put a second, near-identical sequence in the Salmon
    # index that competes with the reference gene for the same reads and dilutes
    # its counts. `o` is left out as ambiguous.
    NOVEL_CLASSES <- c(${params.novel_classes.collect { cc -> "\"${cc}\"" } .join(", ")})
    # isoform_classes are kept in the index but NOT catalogued; their reads are
    # remapped to the reference parent gene downstream (BUILD_TX2GENE).
    ISO_CLASSES   <- c(${params.isoform_classes.collect { cc -> "\"${cc}\"" } .join(", ")})
    if (!("class_code" %in% colnames(rec))) {
        stop("recurrence table lacks class_code — cannot restrict the freeze to novel classes")
    }

    # Filter to high-confidence noncoding that pass recurrence, restricted to
    # genuinely novel classes.
    nc <- subset(rec, consensus_noncoding == TRUE & passes_recurrence == TRUE &
                     !is.na(class_code) & class_code %in% NOVEL_CLASSES)

    if (nrow(nc) == 0) {
        stop("No non-coding transcripts passed all filters. Check your inputs.")
    }

    # Keep these two helpers in sync with BUILD_GENE_CATALOG and BUILD_TX2GENE:
    # both must strip Ensembl versions but leave MSTRG IDs verbatim, or downstream
    # gene_id joins (tx2gene / catalog / DE / report) break silently.
    extract_attr <- function(attrs, key) {
        pat <- paste0('(?:^|; )', key, ' "([^"]+)"')
        m <- regexpr(pat, attrs, perl = TRUE)
        ans <- rep(NA_character_, length(attrs))
        ok <- m > 0L
        if (any(ok)) ans[ok] <- sub(paste0('.*', pat, '.*\$'), '\\\\1', attrs[ok], perl = TRUE)
        ans
    }
    strip_version <- function(x) sub("\\\\.[0-9]+\$", "", x)

    # Frozen transcript ID set (noncoding consensus that pass recurrence)
    nc_ids_requested <- unique(nc\$ID)
    cat("Eligible novel transcripts:", length(nc_ids_requested), "\\n")

    # Extract from the source GTF by EXACT transcript_id attribute. Bare-id
    # `grep -F -f` matches any attribute containing an id as a substring: an
    # eligible reference id stored in cmp_ref drags rejected siblings into the
    # freeze (e.g. MSTRG.5444 matched through cmp_ref ENST00000516707) and
    # prefix collisions (.2 matching .20-.23) leak transcripts that failed
    # recurrence. Anchoring on the quoted transcript_id token selects only
    # lines that belong to eligible transcripts.
    pat <- paste0('transcript_id "', nc_ids_requested, '"')
    writeLines(pat, "frozen_ids.txt")
    system(paste("grep -F -f frozen_ids.txt", shQuote("${source_gtf}"),
                 "> quant_universe.gtf"))

    # Fail loudly instead of silently substituting the full candidate set
    if (file.info("quant_universe.gtf")\$size < 100) {
        stop("No frozen transcripts matched in source GTF. Check transcript ID consistency between CPAT/CPC2 and the assembly.")
    }

    # ---- Drop transcripts without exon structure under their own id ----
    # gffcompare containment ("k"/"c") can store exon features under a sibling
    # transcript_id; the transcript then has no structural fingerprint and would
    # leak catalog_id = NA into DE/report tables. Exclude them here (counted)
    # so frozen, tx2gene and the catalog stay aligned.
    fz_lines <- readLines("quant_universe.gtf")
    fz_parts <- strsplit(fz_lines, "\\t", fixed = TRUE)
    fz_ok    <- lengths(fz_parts) >= 9L
    fz_feat  <- vapply(fz_parts[fz_ok], function(p) p[[3L]], character(1))
    fz_attrs <- vapply(fz_parts[fz_ok], function(p) p[[9L]], character(1))
    fz_txid  <- extract_attr(fz_attrs, "transcript_id")

    frozen_eligible <- unique(na.omit(fz_txid[fz_feat == "transcript"]))
    missing_requested <- setdiff(nc_ids_requested, frozen_eligible)
    if (length(missing_requested) > 0) {
        stop(sprintf("Frozen GTF is missing %d eligible transcript IDs (e.g. %s) — extraction regression",
                     length(missing_requested), paste(head(missing_requested, 5), collapse = ", ")))
    }
    had_exons       <- unique(na.omit(fz_txid[fz_feat == "exon"]))
    no_exons <- setdiff(frozen_eligible, had_exons)
    if (length(no_exons) > 0) {
        cat(sprintf("WARNING: %d eligible transcripts lack exon structure under their own id — excluded from frozen set (e.g. %s)\\n",
                    length(no_exons), paste(head(no_exons, 5), collapse = ", ")))
    }
    drop <- !is.na(fz_txid) & fz_txid %in% no_exons
    fz_lines <- fz_lines[fz_ok][!drop]
    writeLines(fz_lines, "quant_universe.gtf")
    nc_ids <- setdiff(nc_ids_requested, no_exons)
    if (length(nc_ids) == 0) {
        stop("All frozen transcripts lacked exon structure under their own id — check source GTF consistency")
    }

    retained_ids <- unique(na.omit(fz_txid[fz_feat == "transcript" & !drop]))
    if (!setequal(retained_ids, nc_ids)) {
        stop("Frozen GTF does not match the eligible transcript set exactly — extraction regression")
    }

    # ---- Disambiguate IDs reused from the reference ----
    # `stringtie --merge -G` reuses the reference transcript_id for transcripts
    # that only partially match it (class_code != "=", e.g. containment "c"/"k"),
    # which then collides with the same ID in the reference. tximport keeps only
    # the first tx2gene match per transcript, so every read for that transcript
    # silently gets attributed to the reference gene and the novel gene ends up
    # with zero counts everywhere. Disambiguate those IDs before they're used
    # for quantification or tx2gene; exact matches ("=") keep the real accession.
    gtf_lines <- readLines("quant_universe.gtf")
    gtf_ok    <- lengths(strsplit(gtf_lines, "\\t", fixed = TRUE)) >= 9L
    ok_lines  <- gtf_lines[gtf_ok]
    ok_parts  <- strsplit(ok_lines, "\\t", fixed = TRUE)
    feature   <- vapply(ok_parts, function(p) p[[3L]], character(1))
    attrs_all <- vapply(ok_parts, function(p) p[[9L]], character(1))

    tx_mask    <- feature == "transcript"
    tx_ids     <- extract_attr(attrs_all[tx_mask], "transcript_id")
    tx_genes   <- extract_attr(attrs_all[tx_mask], "gene_id")
    tx_class   <- extract_attr(attrs_all[tx_mask], "class_code")

    # Reference transcript id set (version-stripped), for the collision check.
    system(paste("grep -P '\\ttranscript\\t'", shQuote("${reference_gtf}"),
                 "> ref_transcript_lines.txt"))
    ref_lines  <- readLines("ref_transcript_lines.txt")
    ref_parts  <- strsplit(ref_lines, "\\t", fixed = TRUE)
    ref_ok     <- lengths(ref_parts) >= 9L
    ref_attrs  <- vapply(ref_parts[ref_ok], function(p) p[[9L]], character(1))
    ref_id_set <- unique(na.omit(strip_version(extract_attr(ref_attrs, "transcript_id"))))

    needs_rename <- !is.na(tx_class) & tx_class != "=" & strip_version(tx_ids) %in% ref_id_set
    if (any(needs_rename)) {
        old_ids    <- tx_ids[needs_rename]
        dup_old <- unique(old_ids[duplicated(old_ids)])
        if (length(dup_old) > 0) {
            stop(sprintf("Ambiguous disambiguation map: %d reused reference transcript IDs appear under multiple novel entries (e.g. %s)",
                         length(dup_old), paste(head(dup_old, 5), collapse = ", ")))
        }
        rename_map <- setNames(paste0(tx_genes[needs_rename], "__", old_ids), old_ids)

        line_ids <- extract_attr(attrs_all, "transcript_id")
        hit      <- line_ids %in% old_ids
        ok_lines[hit] <- mapply(function(line, old) {
            renamed <- sub(paste0('transcript_id "', old, '"'),
                           paste0('transcript_id "', rename_map[[old]], '"'),
                           line, fixed = TRUE)
            # Keep the original id so downstream consumers (overlap/recurrence
            # evidence joins keyed on the pre-rename id) can re-map to the new id.
            paste0(renamed, sprintf(' orig_transcript_id "%s";', old))
        }, ok_lines[hit], line_ids[hit])

        cat(sprintf("INFO: disambiguated %d transcript IDs reused from the reference (class_code != '=')\\n",
                    length(old_ids)))
        gtf_lines[gtf_ok] <- ok_lines
        writeLines(gtf_lines, "quant_universe.gtf")
    }

    # ---- Split the catalog set from the isoform set ----
    # The quant universe keeps the isoforms (index); the frozen catalog drops
    # them entirely (no transcript/exon lines) so every catalog consumer sees a
    # locus set with no reference-derived isoform in it.
    all_lines <- readLines("quant_universe.gtf")
    parts <- strsplit(all_lines, "\\t", fixed = TRUE)
    okp   <- lengths(parts) >= 9L
    feature <- vapply(parts[okp], function(p) p[[3L]], character(1))
    attrs   <- vapply(parts[okp], function(p) p[[9L]], character(1))
    line_tx <- extract_attr(attrs, "transcript_id")
    line_cl <- extract_attr(attrs, "class_code")
    iso_txids <- unique(na.omit(line_tx[feature == "transcript" &
                                         !is.na(line_cl) & line_cl %in% ISO_CLASSES]))
    is_iso_line <- !is.na(line_tx) & line_tx %in% iso_txids
    writeLines(all_lines[okp][!is_iso_line], "frozen_lncrna.gtf")

    catalog_txids <- setdiff(unique(na.omit(line_tx[feature == "transcript"])), iso_txids)
    n_catalog <- length(catalog_txids)
    n_iso     <- length(iso_txids)
    if (n_catalog == 0) {
        stop("Catalog set is empty after removing isoform_classes — check novel_classes/isoform_classes")
    }
    cat(sprintf("Catalog transcripts: %d; isoforms kept in the index only: %d\\n", n_catalog, n_iso))

    # Checksums
    sha_frozen <- system("sha256sum frozen_lncrna.gtf | cut -d' ' -f1", intern = TRUE)
    sha_quant  <- system("sha256sum quant_universe.gtf | cut -d' ' -f1", intern = TRUE)

    sink("frozen_manifest.json")
    cat(sprintf('{
  "pipeline_version": "0.1.0",
  "freeze_date": "%s",
  "n_novel_lncrna": %d,
  "n_isoform_only": %d,
  "n_quant_universe_transcripts": %d,
  "isoform_classes": [%s],
  "sha256": {
    "frozen_gtf": "%s",
    "quant_universe_gtf": "%s"
  }
}
', Sys.time(), n_catalog, n_iso, n_catalog + n_iso,
         paste(sprintf('"%s"', ISO_CLASSES), collapse = ", "),
         trimws(sha_frozen), trimws(sha_quant)))
    sink()

    writeLines(c(trimws(sha_frozen), trimws(sha_quant)), "frozen_checksums.sha256")
    """
}
