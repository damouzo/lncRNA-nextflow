// Build the novel catalog master (gene-level) and novel transcripts (transcript-level)
// from ASSIGN_CATALOG_IDS + BUILD_GENE_CATALOG + CATALOG_OVERLAP + EXPRESSION_RECURRENCE.
//
// Gene-level join bookkeeping (plan T6, correction 6): after EVERY join the row
// count must match the expected n — we fail fast if any novel gene is lost and
// warn if a join unexpectedly grows the table. Never join by parsing the
// hierarchical catalog ID; always by the parent_catalog_id/catalog_id value
// (correction C).

process BUILD_CATALOG_MASTER {

    input:
    path catalog_ids        // ASSIGN_CATALOG_IDS gene-level (gene_id, catalog_id, locus, class_codes)
    path gene_catalog       // BUILD_GENE_CATALOG full catalog (novel rows carry catalog/ref_gene/cons/synteny)
    path catalog_tx         // ASSIGN_CATALOG_IDS transcript-level map
    path overlap_table      // CATALOG_OVERLAP transcript-level external evidence
    path recurrence_table   // EXPRESSION_RECURRENCE per-transcript metadata
    path frozen_gtf         // frozen novel GTF (transcript<->gene map, exon structure)

    output:
    path "novel_catalog_master.tsv", emit: master
    path "novel_transcripts.tsv",    emit: transcripts

    script:
    """
    #!/usr/bin/env Rscript
    suppressWarnings({
        library(dplyr)
        library(readr)
    })

    extract_attr <- function(attrs, key) {
        pat <- paste0('(?:^|; )', key, ' "([^"]+)"')
        m <- regexpr(pat, attrs, perl = TRUE)
        ans <- rep(NA_character_, length(attrs))
        ok <- m > 0L
        if (any(ok)) ans[ok] <- sub(paste0('.*', pat, '.*\$'), '\\\\1', attrs[ok], perl = TRUE)
        ans
    }

    # ---- Parse frozen GTF (transcript + exon lines) ----
    gtf_lines <- readLines("${frozen_gtf}")
    gtf_parts <- strsplit(gtf_lines, "\\t", fixed = TRUE)
    okp <- lengths(gtf_parts) >= 9L
    parts <- gtf_parts[okp]
    feature <- vapply(parts, function(p) p[[3L]], character(1))
    attrs_all <- vapply(parts, function(p) p[[9L]], character(1))

    f_tx <- data.frame(
        transcript_id = extract_attr(attrs_all[feature == "transcript"], "transcript_id"),
        gene_id       = extract_attr(attrs_all[feature == "transcript"], "gene_id"),
        class_code    = extract_attr(attrs_all[feature == "transcript"], "class_code"),
        orig_transcript_id = extract_attr(attrs_all[feature == "transcript"], "orig_transcript_id"),
        chr    = vapply(parts[feature == "transcript"], `[`, character(1), 1L),
        start  = as.integer(vapply(parts[feature == "transcript"], `[`, character(1), 4L)),
        end    = as.integer(vapply(parts[feature == "transcript"], `[`, character(1), 5L)),
        strand = vapply(parts[feature == "transcript"], `[`, character(1), 7L),
        stringsAsFactors = FALSE)
    f_tx <- f_tx[!is.na(f_tx\$transcript_id) & !is.na(f_tx\$gene_id), ]
    tx2gene <- distinct(f_tx, transcript_id, gene_id)
    frozen_ids <- unique(f_tx\$transcript_id)

    # Transcripts whose ids were disambiguated at freeze carry their original id
    # in `orig_transcript_id`; the overlap/recurrence evidence tables are keyed
    # on that ORIGINAL id, so re-key them to the renamed id before the frozen_ids
    # filter and the tx2gene join keep their evidence.
    has_orig <- !is.na(f_tx\$orig_transcript_id)
    new_to_old <- setNames(f_tx\$orig_transcript_id[has_orig], f_tx\$transcript_id[has_orig])
    remap_tx_id <- function(x) {
        if (length(new_to_old) == 0) return(x)
        old_to_new <- setNames(names(new_to_old), unname(new_to_old))
        y <- as.character(x)
        hit <- y %in% names(old_to_new)
        if (any(hit)) y[hit] <- unname(old_to_new[y[hit]])
        y
    }

    f_ex <- data.frame(
        transcript_id = extract_attr(attrs_all[feature == "exon"], "transcript_id"),
        start  = as.integer(vapply(parts[feature == "exon"], `[`, character(1), 4L)),
        end    = as.integer(vapply(parts[feature == "exon"], `[`, character(1), 5L)),
        stringsAsFactors = FALSE)
    f_ex <- f_ex[!is.na(f_ex\$transcript_id), ]
    ex_struct <- f_ex %>%
        group_by(transcript_id) %>%
        arrange(start) %>%
        summarise(
            n_exons = n(),
            exon_coords = paste(paste0(start, "-", end), collapse = ";"),
            .groups = "drop")

    # ---- Transcript-level external evidence + recurrence ----
    # These tables identify transcripts via the `ID` column (CODING_CONSENSUS
    # naming), NOT `transcript_id`. Renormalize explicitly — inventing an
    # all-NA `transcript_id` column would silently empty the joins below.
    ov <- tryCatch(read_tsv("${overlap_table}", show_col_types = FALSE), error = function(e) data.frame())
    rec <- tryCatch(read_tsv("${recurrence_table}", show_col_types = FALSE), error = function(e) data.frame())
    renorm_tx <- function(df) {
        if (nrow(df) == 0) return(df)
        tcol <- intersect(c("transcript_id", "ID"), colnames(df))
        if (length(tcol) == 0) stop("no transcript ID column in table, cols: ",
                                    paste(colnames(df), collapse = ","))
        if (tcol[1] != "transcript_id") df[["transcript_id"]] <- df[[tcol[1]]]
        df
    }
    ov <- renorm_tx(ov)
    rec <- renorm_tx(rec)
    # Guard evidence/reported columns (always present on this pipeline's CATALOG_OVERLAP
    # output; padding keeps the master stable if a future source ever drops one).
    evidence_cols <- c("known_catalog_hit", "catalog_sources",
                       "lncipedia_id", "noncode_id", "fantomcat_id", "rnacentral_id",
                       "lncipedia_pct", "noncode_pct", "fantomcat_pct", "rnacentral_pct")
    for (cc in setdiff(evidence_cols, colnames(ov))) ov[[cc]] <- NA
    if (nrow(ov) > 0) {
        ov\$transcript_id <- remap_tx_id(ov\$transcript_id)
        ov <- ov[ov\$transcript_id %in% frozen_ids, , drop = FALSE]
    }
    if (nrow(rec) > 0) {
        rec\$transcript_id <- remap_tx_id(rec\$transcript_id)
        rec <- rec[rec\$transcript_id %in% frozen_ids, , drop = FALSE]
    }

    # ---- Gene-level aggregation from transcript tables ----
    cat_tx <- read_tsv("${catalog_tx}", show_col_types = FALSE)

    n_tx_by_gene <- cat_tx %>%
        group_by(parent_catalog_id) %>%
        summarise(n_transcripts = n_distinct(catalog_tx_id), .groups = "drop") %>%
        rename(catalog_id = parent_catalog_id)

    ov_gene <- data.frame(
        gene_id = character(), known_catalog_hit = logical(), catalog_sources = character(),
        lncipedia_id = character(), noncode_id = character(), fantomcat_id = character(),
        rnacentral_id = character(), lncipedia_pct = numeric(), noncode_pct = numeric(),
        fantomcat_pct = numeric(), rnacentral_pct = numeric(), stringsAsFactors = FALSE)
    if (nrow(ov) > 0) {
        ov_gene <- ov %>%
            left_join(tx2gene, by = "transcript_id") %>%
            group_by(gene_id) %>%
            summarise(
                known_catalog_hit = ifelse(any(!is.na(known_catalog_hit) & known_catalog_hit), TRUE, FALSE),
                catalog_sources = paste(unique(na.omit(catalog_sources)), collapse = ","),
                lncipedia_id  = paste(unique(na.omit(lncipedia_id)),  collapse = ","),
                noncode_id    = paste(unique(na.omit(noncode_id)),    collapse = ","),
                fantomcat_id  = paste(unique(na.omit(fantomcat_id)),  collapse = ","),
                rnacentral_id = paste(unique(na.omit(rnacentral_id)), collapse = ","),
                lncipedia_pct  = suppressWarnings(max(0, max(lncipedia_pct,  na.rm = TRUE))),
                noncode_pct    = suppressWarnings(max(0, max(noncode_pct,    na.rm = TRUE))),
                fantomcat_pct  = suppressWarnings(max(0, max(fantomcat_pct,  na.rm = TRUE))),
                rnacentral_pct = suppressWarnings(max(0, max(rnacentral_pct, na.rm = TRUE))),
                .groups = "drop")
        for (cc in c("lncipedia_id", "noncode_id", "fantomcat_id", "rnacentral_id", "catalog_sources")) {
            ov_gene[[cc]][is.na(ov_gene[[cc]]) | ov_gene[[cc]] == ""] <- NA_character_
        }
    }

    rec_gene <- data.frame(
        gene_id = character(), n_reps_expressed = integer(), consensus_noncoding = character(),
        stringsAsFactors = FALSE)
    if (nrow(rec) > 0) {
        rec_gene <- rec %>%
            left_join(tx2gene, by = "transcript_id") %>%
            group_by(gene_id) %>%
            summarise(
                n_reps_expressed = suppressWarnings(max(n_reps_expressed, na.rm = TRUE)),
                consensus_noncoding = "noncoding",
                .groups = "drop") %>%
            mutate(consensus_noncoding = ifelse(is.na(consensus_noncoding), "noncoding", consensus_noncoding))
    }

    assert_nrow <- function(df, expected, label) {
        if (nrow(df) < expected) stop(sprintf("%s: lost rows (have %d, expected %d) — join broken", label, nrow(df), expected))
        if (nrow(df) > expected) cat(sprintf("WARNING: %s grew (have %d, expected %d)\\n", label, nrow(df), expected))
        invisible(df)
    }

    # ---- Gene-level master, with assertion after each join ----
    ids <- read_tsv("${catalog_ids}", show_col_types = FALSE)
    n0 <- nrow(ids)
    if (n0 == 0) { cat("No novel genes to catalog — master empty but valid\\n") }

    gc <- read_tsv("${gene_catalog}", show_col_types = FALSE) %>% filter(origin == "novel")
    if (any(!(ids\$gene_id %in% gc\$gene_id))) stop("BUILD_CATALOG_MASTER: some catalog genes missing from gene_catalog (origin novel)")
    n_genes_no_id <- sum(!(gc\$gene_id %in% ids\$gene_id))
    if (n_genes_no_id > 0) {
        cat(sprintf("WARNING: %d novel genes have no cataloged transcript (all their isoforms lack exon structure) — they are not part of the catalog master\\n",
                    n_genes_no_id))
    }

    m <- left_join(
        select(ids, gene_id, catalog_id, locus, class_codes),
        select(gc, gene_id, gene_name, chromosome, start, end, strand, class_code,
               catalog_id_ = catalog_id, ref_gene_id, ref_gene_name, ref_gene_distance,
               ref_gene_relation, n_exons, conservation_mean, conservation_pct,
               syntenic_locus, syntenic_target_gene_id, u_single_exon_low_confidence),
        by = "gene_id")
    m <- m %>% distinct(gene_id, .keep_all = TRUE)
    assert_nrow(m, n0, "master×gene_catalog")

    m <- left_join(m, n_tx_by_gene, by = "catalog_id")
    assert_nrow(m, n0, "master×n_transcripts")

    m <- left_join(m, ov_gene, by = "gene_id")
    assert_nrow(m, n0, "master×overlap")

    m <- left_join(m, rec_gene, by = "gene_id")
    assert_nrow(m, n0, "master×recurrence")

    m <- m %>%
        mutate(
            class_codes = coalesce(class_codes, class_code),
            across(c(lncipedia_id, noncode_id, fantomcat_id, rnacentral_id),
                   ~ ifelse(is.na(.x) | !nzchar(.x), NA_character_, .x)),
            catalog_sources = ifelse(is.na(catalog_sources) | !nzchar(catalog_sources), NA_character_, catalog_sources),
            known_catalog_hit = ifelse(is.na(known_catalog_hit), FALSE, known_catalog_hit)) %>%
        select(catalog_id, gene_id, gene_name, locus, chromosome, start, end, strand,
               n_transcripts, n_exons, class_codes, consensus_noncoding, n_reps_expressed,
               ref_gene_id, ref_gene_name, ref_gene_distance, ref_gene_relation,
               known_catalog_hit, catalog_sources,
               lncipedia_id, lncipedia_pct, noncode_id, noncode_pct,
               fantomcat_id, fantomcat_pct, rnacentral_id, rnacentral_pct,
               conservation_mean, conservation_pct, syntenic_locus, syntenic_target_gene_id,
               u_single_exon_low_confidence)
    assert_nrow(m, n0, "master-final")

    if (any(is.na(m\$catalog_id))) stop("catalog_id NA for a novel gene in master")

    write_tsv(m, "novel_catalog_master.tsv")

    # ---- Transcript-level table (exon/intron structure) ----
    tr <- cat_tx %>%
        mutate(transcript_id = mstrg_tx_id) %>%
        left_join(tx2gene, by = "transcript_id") %>%
        left_join(select(f_tx, transcript_id, chr, start, end, strand, class_code), by = "transcript_id") %>%
        left_join(ex_struct, by = "transcript_id") %>%
        left_join(select(ov, transcript_id, known_catalog_hit,
                         lncipedia_id, noncode_id, fantomcat_id, rnacentral_id), by = "transcript_id") %>%
        mutate(known_catalog_hit = ifelse(is.na(known_catalog_hit), FALSE, known_catalog_hit))
    if (any(is.na(tr\$transcript_id))) stop("catalog transcript not found in the frozen GTF")
    assert_nrow(tr, nrow(cat_tx), "transcripts")
    # The catalog universe is defined by ASSIGN_CATALOG_IDS (catalog_tx), not by
    # the raw frozen GTF: transcripts whose exon features carry a different
    # transcript_id (gffcompare "k"/"c" containment) have no structural
    # signature and are intentionally absent from the registry. Warn, never fail.
    n_skipped <- length(frozen_ids) - nrow(cat_tx)
    if (n_skipped > 0) {
        cat(sprintf("WARNING: %d frozen transcripts lack exon structure under their own id — excluded from the catalog (see ASSIGN_CATALOG_IDS)\\n",
                    n_skipped))
    }

    write_tsv(tr, "novel_transcripts.tsv")
    cat("Master:", nrow(m), "genes; transcripts:", nrow(tr), "\\n")
    """
}