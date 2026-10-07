// Gene catalog of the full analysis universe
// One row per gene_id, tagging each as:
//   novel            — frozen discovery lncRNAs (with gffcompare class_code)
//   annotated_lncrna — reference genes whose gene_biotype is in the configurable list
//   reference_other  — everything else in the reference annotation
// Used downstream (reporting, cis, enrichment) without touching quantification.

process BUILD_GENE_CATALOG {

    input:
    path frozen_gtf         // frozen novel lncRNA GTF (class_code in attributes)
    path reference_gtf      // reference annotation
    path conservation_scores  // per-transcript conservation (header-only placeholder when skipped)
    path synteny_scores       // per-transcript synteny (header-only placeholder when skipped)
    path catalog_ids          // ASSIGN_CATALOG_IDS gene-level table (catalog_id, locus)
    path overlap_table        // CATALOG_OVERLAP transcript-level external hits (reporting)

    output:
    path "gene_catalog.tsv", emit: gene_catalog
    path "catalog_stats.tsv", emit: catalog_stats

    script:
    """
    #!/usr/bin/env Rscript
    library(dplyr)
    library(readr)

    ANNOT_BIOTYPES <- c(${params.annotated_lncrna_biotypes.collect { bt -> "\"${bt}\"" } .join(", ")})

    extract_attr <- function(attrs, key) {
        pat <- paste0('(?:^|; )', key, ' "([^"]+)"')
        m <- regexpr(pat, attrs, perl = TRUE)
        ans <- rep(NA_character_, length(attrs))
        ok <- m > 0L
        if (any(ok)) ans[ok] <- sub(paste0('.*', pat, '.*\$'), '\\\\1', attrs[ok], perl = TRUE)
        ans
    }
    strip_version <- function(x) sub("\\\\.[0-9]+\$", "", x)

    fill_table <- function(lines) {
        if (length(lines) == 0) {
            return(data.frame(
                chromosome = character(),
                start      = integer(),
                end        = integer(),
                strand     = character(),
                attrs      = character(),
                stringsAsFactors = FALSE
            ))
        }
        parts <- strsplit(lines, "\\t", fixed = TRUE)
        idx <- lengths(parts) >= 9L
        if (!any(idx)) stop("No parseable lines in this GTF subset")
        parts <- parts[idx]
        data.frame(
            chromosome = vapply(parts, `[`, character(1), 1L),
            start      = as.integer(vapply(parts, `[`, character(1), 4L)),
            end        = as.integer(vapply(parts, `[`, character(1), 5L)),
            strand     = vapply(parts, `[`, character(1), 7L),
            attrs      = vapply(parts, `[`, character(1), 9L),
            stringsAsFactors = FALSE
        )
    }

    # ---- Reference genes (one line per gene) ----
    system(paste("grep -P '\\tgene\\t'", shQuote("${reference_gtf}"),
                 "> ref_gene_lines.txt"))
    ref <- fill_table(readLines("ref_gene_lines.txt"))
    ref\$gene_id      <- strip_version(extract_attr(ref\$attrs, "gene_id"))
    ref\$gene_name    <- extract_attr(ref\$attrs, "gene_name")
    ref\$gene_biotype <- extract_attr(ref\$attrs, "gene_biotype")

    ref_cat <- ref %>%
        filter(!is.na(gene_id)) %>%
        mutate(
            origin     = ifelse(gene_biotype %in% ANNOT_BIOTYPES, "annotated_lncrna", "reference_other"),
            class_code = NA_character_,
            # GRanges only accepts + - *; some GTFs use "." or empty for unknown strand
            strand     = ifelse(strand %in% c("+", "-", "*"), strand, "*")
        ) %>%
        distinct(gene_id, .keep_all = TRUE) %>%
        select(gene_id, gene_name, chromosome, start, end, strand,
               gene_biotype, origin, class_code)

    # ---- Novel lncRNAs (transcript-level, aggregated to gene) ----
    # MSTRG IDs are not versioned like Ensembl — keep them verbatim
    system(paste("grep -P '\\ttranscript\\t'", shQuote("${frozen_gtf}"),
                 "> novel_transcript_lines.txt"))
    nov_lines <- readLines("novel_transcript_lines.txt")
    if (length(nov_lines) == 0) {
        cat("WARNING: no novel transcripts parsed from frozen GTF, novel origin will be empty\\n")
    }
    nov <- fill_table(nov_lines)
    nov\$gene_id <- extract_attr(nov\$attrs, "gene_id")
    nov\$class_code <- extract_attr(nov\$attrs, "class_code")

    nov_cat <- nov %>%
        filter(!is.na(gene_id)) %>%
        group_by(gene_id) %>%
        summarise(
            chromosome = dplyr::first(chromosome),
            start      = min(start),
            end        = max(end),
            strand     = ifelse(first(strand) %in% c("+", "-", "*"), first(strand), "*"),
            class_code = paste(sort(unique(na.omit(class_code))), collapse = ","),
            .groups    = "drop"
        ) %>%
        mutate(
            origin       = "novel",
            gene_biotype = NA_character_,
            gene_name    = NA_character_
        ) %>%
        select(gene_id, gene_name, chromosome, start, end, strand,
               gene_biotype, origin, class_code)

    # ---- Merge: novel wins if a gene_id ever collides with reference ----
    catalog <- bind_rows(nov_cat, ref_cat) %>%
        distinct(gene_id, .keep_all = TRUE) %>%
        arrange(origin, gene_id)

    # ---- Conservation & synteny (reporting-only; all-NA when skipped) ----
    # Scores are transcript-level for the frozen (novel) set. Aggregate to the
    # gene level using the transcript->gene map from the frozen GTF. Genes
    # without data (reference, or modules skipped) stay NA — never 0, so the
    # report can tell "not evaluated" apart from "not conserved".
    map_tx2gene <- function(lines) {
        if (length(lines) == 0) return(data.frame(transcript_id = character(),
                                                  gene_id = character()))
        parts <- strsplit(lines, "\\t", fixed = TRUE)
        attrs <- vapply(parts, function(p) p[[9L]], character(1))
        data.frame(
            transcript_id = extract_attr(attrs, "transcript_id"),
            gene_id       = extract_attr(attrs, "gene_id"),
            orig_transcript_id = extract_attr(attrs, "orig_transcript_id"),
            stringsAsFactors = FALSE
        )
    }
    system(paste("grep -P '\\ttranscript\\t'", shQuote("${frozen_gtf}"),
                 "> tx2gene_frozen_lines.txt"))
    tx2gene_frozen <- map_tx2gene(readLines("tx2gene_frozen_lines.txt")) %>%
        filter(!is.na(transcript_id), !is.na(gene_id)) %>%
        distinct(transcript_id, .keep_all = TRUE)

    # Transcripts whose ids were disambiguated at freeze (reference id reused with
    # class_code != "=") carry the original id in `orig_transcript_id`; external
    # catalog evidence (CATALOG_OVERLAP) is keyed on that ORIGINAL id and must be
    # re-keyed to the renamed id before the join below.
    has_orig <- !is.na(tx2gene_frozen\$orig_transcript_id)
    new_to_old <- setNames(tx2gene_frozen\$orig_transcript_id[has_orig],
                           tx2gene_frozen\$transcript_id[has_orig])
    remap_tx_id <- function(x) {
        if (length(new_to_old) == 0) return(x)
        old_to_new <- setNames(names(new_to_old), unname(new_to_old))
        y <- as.character(x)
        hit <- y %in% names(old_to_new)
        if (any(hit)) y[hit] <- unname(old_to_new[y[hit]])
        y
    }

read_scores <- function(f, cols) {
    x <- tryCatch(read_tsv(f, show_col_types = FALSE), error = function(e) NULL)
    if (is.null(x)) {
        # Unreadable file: build a 0-row frame with the expected columns so
        # downstream joins/summaries still resolve their symbols.
        x <- as.data.frame(setNames(lapply(cols, function(cc) character(0)), cols),
                           stringsAsFactors = FALSE)
    } else {
        x <- x[x\$transcript_id %in% tx2gene_frozen\$transcript_id, ]
    }
    # For a header-only placeholder (0 rows) the columns are already correct;
    # they just aggregate to NA downstream.
    x
}

cons  <- read_scores("${conservation_scores}",
                     c("transcript_id", "mean_score", "max_score", "pct_bases_conserved"))
synth <- read_scores("${synteny_scores}",
                     c("transcript_id", "syntenic_locus", "syntenic_target_gene_id"))

    catalog <- catalog %>%
        left_join(
            cons %>%
                left_join(tx2gene_frozen, by = "transcript_id") %>%
                group_by(gene_id) %>%
                summarise(
                    conservation_mean = suppressWarnings(mean(as.numeric(mean_score), na.rm = TRUE)),
                    conservation_pct  = suppressWarnings(mean(as.numeric(pct_bases_conserved), na.rm = TRUE)),
                    .groups = "drop"
                ),
            by = "gene_id"
        ) %>%
        left_join(
            synth %>%
                left_join(tx2gene_frozen, by = "transcript_id") %>%
                group_by(gene_id) %>%
                summarise(
                    syntenic_locus = ifelse(any(tolower(syntenic_locus) == "true"), TRUE,
                                           ifelse(all(is.na(syntenic_locus)), NA, FALSE)),
                    syntenic_target_gene_id = paste(unique(na.omit(syntenic_target_gene_id)),
                                                    collapse = ","),
                    .groups = "drop"
                ) %>%
                mutate(syntenic_target_gene_id = ifelse(nchar(syntenic_target_gene_id) == 0,
                                                        NA_character_, syntenic_target_gene_id)),
            by = "gene_id"
        )

    # ---- Monoexonic u (lincRNA) low-confidence flag ----
    # Count exons per gene_id in the frozen GTF, then flag novel genes with
    # class_code "u" and a single exon. Reporting-only — never a filter.
    system(paste("grep -P '\\texon\\t'", shQuote("${frozen_gtf}"),
                 "> frozen_exon_lines.txt"))
    exon_lines <- readLines("frozen_exon_lines.txt")
    if (length(exon_lines) > 0) {
        exon <- fill_table(exon_lines)
        exon\$gene_id <- extract_attr(exon\$attrs, "gene_id")
        n_exons_per_gene <- exon %>%
            filter(!is.na(gene_id)) %>%
            distinct(gene_id, start, end) %>%
            count(gene_id, name = "n_exons")
    } else {
        n_exons_per_gene <- data.frame(gene_id = character(), n_exons = integer())
    }

    catalog <- catalog %>%
        left_join(n_exons_per_gene, by = "gene_id") %>%
        mutate(
            n_exons = coalesce(n_exons, 0L),
            u_single_exon_low_confidence = origin == "novel" &
                grepl("u", class_code, fixed = TRUE) &
                n_exons == 1L
        )

    # ---- Stable catalog identity + associated reference gene ----
    # catalog_ids.tsv is produced by ASSIGN_CATALOG_IDS (gene level). It may be
    # absent on runs with the registry feature disabled — degrade to NA silently.
    cat_ids <- tryCatch(read_tsv("${catalog_ids}", show_col_types = FALSE), error = function(e) NULL)
    if (!is.null(cat_ids) && "gene_id" %in% colnames(cat_ids) &&
        all(c("catalog_id", "locus") %in% colnames(cat_ids))) {
        # Do NOT pre-initialize catalog_id/locus on `catalog`: both sides would
        # carry them and left_join() would suffix the shared columns (.x/.y),
        # dropping the plain catalog_id the master later selects.
        catalog <- left_join(catalog, select(cat_ids, gene_id, catalog_id, locus), by = "gene_id")
    } else {
        catalog\$catalog_id <- NA_character_
        catalog\$locus      <- NA_character_
    }

    # ---- Reference protein-coding genes: nearest / overlapping / host ----
    refcoding <- ref %>%
        filter(!is.na(gene_id), gene_biotype == "protein_coding",
               !is.na(chromosome), !is.na(start), !is.na(end)) %>%
        distinct(gene_id, gene_name, chromosome, start, end, strand)

    which.max <- function(x) which(x == max(x, na.rm = TRUE))[1]
    which.min <- function(x) which(x == min(x, na.rm = TRUE))[1]

    assign_ref <- function(chr, start, end, strand, refcoding) {
        base <- list(id = NA_character_, name = NA_character_, dist = NA_integer_, rel = NA_character_)
        if (is.na(chr) || is.na(start) || is.na(end)) return(base)
        sub <- refcoding[refcoding\$chromosome == chr, , drop = FALSE]
        if (nrow(sub) == 0) return(list(id = base\$id, name = base\$name, dist = base\$dist, rel = "intergenic"))
        w <- pmax(0L, pmin(end, sub\$end) - pmax(start, sub\$start) + 1L)
        ov <- which(w > 0L)
        if (length(ov) > 0) {
            best <- ov[which.max(w[ov])]
            same <- substr(sub\$strand[best], 1, 1) == strsplit(strand, "")[[1]][1]
            contained <- start >= sub\$start[best] && end <= sub\$end[best]
            rel <- if (same) { if (contained) "intronic" else "overlap_sense" } else "overlap_antisense"
            return(list(id = sub\$gene_id[best], name = sub\$gene_name[best], dist = 0L, rel = rel))
        }
        dist <- pmin(abs(start - sub\$end), abs(sub\$start - end))
        k <- which.min(dist)
        rel <- if (end <= sub\$start[k]) "upstream" else if (start >= sub\$end[k]) "downstream" else "intergenic"
        list(id = sub\$gene_id[k], name = sub\$gene_name[k], dist = as.integer(dist[k]), rel = rel)
    }

    catalog\$ref_gene_id       <- NA_character_
    catalog\$ref_gene_name     <- NA_character_
    catalog\$ref_gene_distance <- NA_integer_
    catalog\$ref_gene_relation <- NA_character_
    nov_idx <- which(catalog\$origin == "novel")
    for (i in nov_idx) {
        rr <- assign_ref(catalog\$chromosome[i], catalog\$start[i], catalog\$end[i],
                         catalog\$strand[i], refcoding)
        catalog\$ref_gene_id[i]       <- rr\$id
        catalog\$ref_gene_name[i]     <- rr\$name
        catalog\$ref_gene_distance[i] <- rr\$dist
        catalog\$ref_gene_relation[i] <- rr\$rel
    }

    # ---- External catalog hits (reporting): aggregate transcript-level evidence ----
    ov_tab <- tryCatch(read_tsv("${overlap_table}", show_col_types = FALSE), error = function(e) data.frame())
    # CATALOG_OVERLAP keeps the CODING_CONSENSUS id column (`ID`), not
    # `transcript_id` — renormalize like build_catalog_master's renorm_tx().
    tid_over <- intersect(c("transcript_id", "ID"), colnames(ov_tab))
    if (nrow(ov_tab) > 0 && length(tid_over) == 1) ov_tab\$transcript_id <- ov_tab[[tid_over[1]]]
    if (nrow(ov_tab) > 0 && "transcript_id" %in% colnames(ov_tab) &&
        "known_catalog_hit" %in% colnames(ov_tab)) {
        ov_tab\$transcript_id <- remap_tx_id(ov_tab\$transcript_id)
        ov_keep <- ov_tab[ov_tab\$transcript_id %in% tx2gene_frozen\$transcript_id, , drop = FALSE]
        ov_gene <- ov_keep %>%
            left_join(tx2gene_frozen, by = "transcript_id") %>%
            group_by(gene_id) %>%
            summarise(
                known_catalog_hit = any(!is.na(known_catalog_hit) & known_catalog_hit),
                catalog_sources = paste(unique(na.omit(catalog_sources)), collapse = ","),
                .groups = "drop")
        ov_gene\$known_catalog_hit[is.na(ov_gene\$known_catalog_hit)] <- FALSE
        catalog <- left_join(catalog, ov_gene, by = "gene_id") %>%
            mutate(
                known_catalog_hit = ifelse(origin == "novel", coalesce(known_catalog_hit, FALSE), NA),
                catalog_sources   = ifelse(origin == "novel" & !is.na(catalog_sources), catalog_sources, NA_character_)
            )
    } else {
        catalog\$known_catalog_hit <- NA
        catalog\$catalog_sources   <- NA_character_
    }

    if (nrow(catalog) == 0) stop("Gene catalog is empty — check reference GTF parsing")

    write_tsv(catalog, "gene_catalog.tsv")

    stats <- catalog %>%
        count(origin, name = "n_genes") %>%
        bind_rows(
            catalog %>% filter(!is.na(gene_biotype)) %>%
                count(gene_biotype, name = "n_genes") %>%
                rename(origin = gene_biotype)
        )
    write_tsv(stats, "catalog_stats.tsv")

    cat("Catalog genes:", nrow(catalog), "\\n")
    print(table(catalog\$origin, useNA = "ifany"))
    """
}