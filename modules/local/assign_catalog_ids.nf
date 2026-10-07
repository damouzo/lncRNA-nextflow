// Assign stable hierarchical catalog IDs (gene + transcript) backed by a
// persistent two-level registry outside results/.
//
// The registry is a TSV with a `entity_type` column ('gene'|'transcript').
// It is never regenerated: it is read and extended in place, atomically.
// The catalog_id of a novel gene is resolved through its transcripts first
// (inherited from matched transcripts via parent_catalog_id); only when no
// transcript matches is the gene-level span fallback evaluated. This ordering
// is what keeps the catalog_id stable when a novel isoform extends the gene
// span (the pre-existing transcripts still match by structural signature).
//
// `cache false` is a deliberate exception: this process writes outside the
// Nextflow work directory (the registry), which breaks the -resume idempotency
// contract. Do not "optimize" it back on.

process ASSIGN_CATALOG_IDS {

    label 'process_medium'
    memory 24.GB
    cache false

    input:
    path frozen_gtf
    val  catalog_registry            // path may not exist yet — passed as val, not file()
    val  genome_build                // fail-fast guard: registry annotated for a different build
    val  reference_annotation        // GTF path; only the basename is recorded
    val  catalog_registry_project    // provenance label

    output:
    path "frozen_lncrna.catalog.gtf", emit: catalog_gtf
    path "catalog_ids.tsv",           emit: catalog_ids
    path "catalog_tx.tsv",            emit: catalog_tx
    path "registry_review_needed.tsv", optional: true, emit: registry_review

    script:
    """
    #!/usr/bin/env Rscript
    suppressWarnings({
        library(dplyr)
        library(readr)
    })

    PREFIX   <- "${params.catalog_id_prefix}"
    MIN_REC  <- as.numeric("${params.catalog_registry_min_reciprocal}")
    BUILD    <- "${genome_build}"
    if (is.null(BUILD) || !nzchar(BUILD) || BUILD %in% c("null", "[]", "none")) BUILD <- ""
    REFANN   <- basename("${reference_annotation}")
    PROJ     <- "${catalog_registry_project}"
    if (is.null(PROJ) || !nzchar(PROJ) || PROJ %in% c("null", "[]")) PROJ <- ""
    RUN_TOKEN <- format(Sys.time(), "%Y%m%dT%H%M%S")
    nn <- function(x) !is.null(x) && length(x) > 0 && !is.na(x) && nzchar(x)

    reg_path <- "${catalog_registry}"
    reg_off  <- is.null(reg_path) || is.na(reg_path) || !nzchar(reg_path) ||
                reg_path %in% c("null", "[]", "none")
    if (!reg_off && !grepl("^/", reg_path)) {
        # The registry lives outside the Nextflow work dir; a relative path would
        # silently resolve against this task's scratch dir and break persistence.
        stop("catalog_registry must be an ABSOLUTE path, got: ", reg_path)
    }

    # ---- GTF attribute helpers (mirror BUILD_GENE_CATALOG) ----
    extract_attr <- function(attrs, key) {
        pat <- paste0('(?:^|; )', key, ' "([^"]+)"')
        m <- regexpr(pat, attrs, perl = TRUE)
        ans <- rep(NA_character_, length(attrs))
        ok <- m > 0L
        if (any(ok)) ans[ok] <- sub(paste0('.*', pat, '.*\$'), '\\\\1', attrs[ok], perl = TRUE)
        ans
    }

    # ---- Parse frozen GTF ----
    gtf_lines <- readLines("${frozen_gtf}")
    gtf_parts <- strsplit(gtf_lines, "\\t", fixed = TRUE)
    okp <- lengths(gtf_parts) >= 9L
    parts <- gtf_parts[okp]
    feature   <- vapply(parts, function(p) p[[3L]], character(1))
    attrs_all <- vapply(parts, function(p) p[[9L]], character(1))

    tx_ok <- feature == "transcript"
    tx_dat <- data.frame(
        transcript_id = extract_attr(attrs_all[tx_ok], "transcript_id"),
        gene_id       = extract_attr(attrs_all[tx_ok], "gene_id"),
        class_code    = extract_attr(attrs_all[tx_ok], "class_code"),
        chr    = vapply(parts[tx_ok], `[`, character(1), 1L),
        start  = as.integer(vapply(parts[tx_ok], `[`, character(1), 4L)),
        end    = as.integer(vapply(parts[tx_ok], `[`, character(1), 5L)),
        strand = vapply(parts[tx_ok], `[`, character(1), 7L),
        stringsAsFactors = FALSE)
    tx_dat <- tx_dat[!is.na(tx_dat\$transcript_id) & !is.na(tx_dat\$gene_id), ]

    exon_ok <- feature == "exon"
    ex_dat <- data.frame(
        transcript_id = extract_attr(attrs_all[exon_ok], "transcript_id"),
        chr    = vapply(parts[exon_ok], `[`, character(1), 1L),
        start  = as.integer(vapply(parts[exon_ok], `[`, character(1), 4L)),
        end    = as.integer(vapply(parts[exon_ok], `[`, character(1), 5L)),
        stringsAsFactors = FALSE)
    ex_dat <- ex_dat[!is.na(ex_dat\$transcript_id), ]

    if (nrow(tx_dat) == 0) stop("No transcript features parsed from frozen GTF")

    # ---- Structural signature per transcript ----
    # Multiexonic: ordered chain of intron coordinates + strand.
    # Monoexonic: exon extent + strand.
    sig_v   <- character(nrow(tx_dat))
    coords_v <- character(nrow(tx_dat))
    span_s  <- integer(nrow(tx_dat))
    span_e  <- integer(nrow(tx_dat))
    n_ex    <- integer(nrow(tx_dat))
    for (i in seq_len(nrow(tx_dat))) {
        ex <- ex_dat[ex_dat\$transcript_id == tx_dat\$transcript_id[i], ]
        ex <- ex[order(ex\$start), , drop = FALSE]
        n_ex[i] <- nrow(ex)
        if (nrow(ex) == 0) {
            sig_v[i] <- NA_character_; coords_v[i] <- NA_character_
            span_s[i] <- NA_integer_; span_e[i] <- NA_integer_
        } else if (nrow(ex) == 1) {
            sig_v[i] <- paste0(tx_dat\$chr[i], "|", ex\$start[1], "-", ex\$end[1], "|", tx_dat\$strand[i])
            coords_v[i] <- paste0(ex\$start[1], "-", ex\$end[1])
            span_s[i] <- ex\$start[1]; span_e[i] <- ex\$end[1]
        } else {
            int_d <- paste(mapply(function(a, b) paste0(a, "-", b),
                                  ex\$end[-nrow(ex)] + 1L, ex\$start[-1L] - 1L), collapse = ";")
            sig_v[i] <- paste0(tx_dat\$chr[i], "|", int_d, "|", tx_dat\$strand[i])
            coords_v[i] <- paste(paste0(ex\$start, "-", ex\$end), collapse = ";")
            span_s[i] <- min(ex\$start); span_e[i] <- max(ex\$end)
        }
    }
    tx_dat\$signature   <- sig_v
    tx_dat\$exon_coords <- coords_v
    tx_dat\$span_start  <- span_s
    tx_dat\$span_end    <- span_e
    tx_dat\$n_exons     <- n_ex
    n_no_sig <- nrow(tx_dat) - sum(!is.na(tx_dat\$signature))
    tx_dat <- tx_dat[!is.na(tx_dat\$signature), ]
    if (n_no_sig > 0) {
        # gffcompare "k"/"c" transcripts can have their exon features under a
        # sibling/containing transcript_id — they are not registered (no stable
        # structural fingerprint), but stay part of the frozen run.
        cat(sprintf("WARNING: %d frozen transcripts have no exon structure under their own id — excluded from the registry/catalog\\n",
                    n_no_sig))
    }

    # ---- Registry read / parse ----
    REG_COLS <- c("entity_type", "catalog_id", "parent_catalog_id", "mstrg_id",
                  "chr", "start", "end", "strand", "class_code",
                  "splice_junction_signature", "exon_coords",
                  "first_seen_project", "first_seen_run", "first_seen_date")

    read_registry <- function(path) {
        lines <- readLines(path)
        meta  <- list()
        comments <- lines[grepl("^#", lines)]
        for (c in comments) {
            # Regex (not fixed=TRUE, which treated "^#" as literal): strips any
            # leading "#"/space mix so meta keys parse and the fail-fast guards work.
            cc <- sub("^[#[:space:]]+", "", c)
            if (grepl("=", cc, fixed = TRUE)) {
                kv <- strsplit(cc, "=", fixed = TRUE)[[1]]
                if (length(kv) >= 2) meta[[trimws(kv[1])]] <- trimws(paste(kv[-1], collapse = "="))
            }
        }
        dat <- lines[!grepl("^#", lines) & nzchar(lines)]
        if (length(dat) == 0) return(list(meta = meta, data = NULL))
        hdr  <- trimws(strsplit(dat[1], "\\t", fixed = TRUE)[[1]])
        body <- strsplit(dat[-1], "\\t", fixed = TRUE)
        body <- body[lengths(body) >= length(hdr)]
        if (length(body) == 0) return(list(meta = meta, data = NULL))
        m <- matrix(NA_character_, nrow = length(body), ncol = length(hdr))
        for (i in seq_along(body)) m[i, seq_len(min(length(body[[i]]), length(hdr)))] <- body[[i]][seq_len(min(length(body[[i]]), length(hdr)))]
        df <- as.data.frame(m, stringsAsFactors = FALSE)
        colnames(df) <- hdr
        for (cc in colnames(df)) {
            if (is.character(df[[cc]])) {
                df[[cc]] <- trimws(df[[cc]])
                df[[cc]][is.na(df[[cc]]) | df[[cc]] == "NA"] <- NA_character_
            }
        }
        if ("start" %in% colnames(df)) df\$start <- suppressWarnings(as.numeric(df\$start))
        if ("end"   %in% colnames(df)) df\$end   <- suppressWarnings(as.numeric(df\$end))
        list(meta = meta, data = df)
    }

    write_registry <- function(tmp, meta, df) {
        # `tmp` is the FINAL temporary path (no ".tmp" appended here — the
        # caller already adds it). Writing to this exact path is what lets the
        # atomic file.rename() below succeed.
        if (file.exists(tmp)) unlink(tmp)
        sink(tmp)
        # Close the redirection even on error inside this function.
        on.exit(if (sink.number() > 0) sink(), add = TRUE)
        for (nm in names(meta)) cat(sprintf("# %s=%s\\n", nm, meta[[nm]]))
        # Explicit sep = "": the default " " would leave a trailing space on line.
        cat(paste(REG_COLS, collapse = "\\t"), "\\n", sep = "")
        if (!is.null(df) && nrow(df) > 0) {
            missing <- setdiff(REG_COLS, colnames(df))
            if (length(missing) > 0) stop("registry missing columns: ", paste(missing, collapse = ","))
            for (i in seq_len(nrow(df))) {
                cat(paste(as.character(df[i, REG_COLS]), collapse = "\\t"), "\\n", sep = "")
            }
        }
        sink()
    }

    # ---- Schema alignment: trim names/values, reorder, pad missing REG_COLS ----
    normalize_cols <- function(df) {
        if (is.null(df)) return(df)
        if (nrow(df) == 0) {
            out <- data.frame(matrix(NA_character_, nrow = 0, ncol = length(REG_COLS)))
            colnames(out) <- REG_COLS
            return(out)
        }
        colnames(df) <- trimws(colnames(df))
        for (cc in intersect(REG_COLS, colnames(df)))
            if (is.character(df[[cc]])) df[[cc]] <- trimws(df[[cc]])
        for (cc in setdiff(REG_COLS, colnames(df))) df[[cc]] <- NA_character_
        df[, REG_COLS]
    }

    # ---- Interval overlap helpers ----
    parse_exons <- function(s) {
        if (is.na(s) || !nzchar(s)) return(matrix(0L, nrow = 0, ncol = 2))
        toks <- strsplit(s, ";", fixed = TRUE)[[1]]
        toks <- toks[nzchar(toks)]
        if (length(toks) == 0) return(matrix(0L, nrow = 0, ncol = 2))
        do.call(rbind, lapply(toks, function(t) {
            p <- as.integer(strsplit(t, "-", fixed = TRUE)[[1]])
            as.integer(p[1:2])
        }))
    }
    intervals_isec <- function(A, B) {
        if (nrow(A) == 0 || nrow(B) == 0) return(0L)
        tot <- 0L
        for (i in seq_len(nrow(A))) for (j in seq_len(nrow(B))) {
            lo <- max(A[i, 1], B[j, 1]); hi <- min(A[i, 2], B[j, 2])
            if (hi >= lo) tot <- tot + (hi - lo + 1L)
        }
        tot
    }
    exonic_recip <- function(qcoords, ccoords) {
        A <- parse_exons(qcoords); B <- parse_exons(ccoords)
        isec <- intervals_isec(A, B)
        if (isec == 0) return(0)
        qlen <- sum(A[, 2] - A[, 1] + 1L); clen <- sum(B[, 2] - B[, 1] + 1L)
        min(isec / qlen, isec / clen)
    }
    span_recip <- function(s1, e1, s2, e2) {
        lo <- max(s1, s2); hi <- min(e1, e2)
        ov <- if (hi >= lo) (hi - lo + 1L) else 0L
        if (ov == 0) return(0)
        w1 <- e1 - s1 + 1L
        w2 <- e2 - s2 + 1L
        min(ov / w1, ov / w2)
    }

    gene_num <- function(x) {
        if (is.na(x) || !nzchar(x)) return(0L)
        b <- sub(PREFIX, "", x, fixed = TRUE)
        suppressWarnings(as.integer(strsplit(b, ".", fixed = TRUE)[[1]][1]))
    }
    tx_suffix <- function(x) {
        if (is.na(x) || !nzchar(x)) return(NA_integer_)
        b <- sub(PREFIX, "", x, fixed = TRUE)
        p <- strsplit(b, ".", fixed = TRUE)[[1]]
        suppressWarnings(as.integer(p[length(p)]))
    }

    # ---- Main registry write: lock + resolve + persist ----
    # Registered in a function frame so on.exit() fires reliably on error
    # (a top-level on.exit in a script does not).
    run_registry_write <- function() {
        reg_gene <- data.frame()
        reg_tx   <- data.frame()
        reg_df   <- NULL
        meta     <- list()
        lock     <- ""
        global_max <- 0L

        if (reg_off) {
            cat("WARNING: catalog_registry is null — catalog IDs are ephemeral (feature off).\\n")
        } else {
            rdir <- dirname(reg_path)
            if (!nzchar(rdir) || rdir == ".") rdir <- getwd()
            if (!dir.exists(rdir)) dir.create(rdir, recursive = TRUE)
            if (!dir.exists(file.path(rdir, "history"))) dir.create(file.path(rdir, "history"), recursive = TRUE)

            # Atomic lock via mkdir (no pre-check: fails only if it already exists ->
            # TOCTOU-safe). An `owner` marker (PID + node + timestamp) distinguishes a
            # live concurrent run from an orphaned one; stale/orphaned locks are
            # removed automatically so an interrupted run (e.g. OOM-kill before the
            # on.exit could run) can never wedge future runs. The lock is acquired
            # BEFORE the registry is read so the read-modify-write cycle is
            # serialized: a run that starts while another one holds the lock can
            # never persist rows computed from a pre-lock snapshot (which would
            # silently drop the earlier run's committed rows and churn catalog ids).
            LOCK_STALE_SECS <- 1800
            lock <- paste0(reg_path, ".lock")
            owner_file <- file.path(lock, "owner")
            acquired <- FALSE
            for (attempt in 1:2) {
                if (dir.create(lock, recursive = FALSE)) { acquired <- TRUE; break }
                stale <- FALSE
                if (!file.exists(owner_file)) {
                    stale <- TRUE
                } else {
                    mtime <- suppressWarnings(as.numeric(file.info(owner_file)\$mtime))
                    if (is.na(mtime) || (as.numeric(Sys.time()) - mtime) > LOCK_STALE_SECS) stale <- TRUE
                }
                if (!stale) break
                cat("Removing stale catalog registry lock (orphaned): ", lock, "\\n")
                unlink(lock, recursive = TRUE)
            }
            if (!acquired) {
                owner <- if (file.exists(owner_file)) paste(readLines(owner_file, warn = FALSE), collapse = " ") else "unknown"
                stop("catalog registry lock exists and looks live: ", lock,
                     " (owner: ", owner, "). A concurrent run may be writing — wait for it, ",
                     "or remove the lock manually if it is truly orphaned.")
            }
            # Function-frame on.exit: runs on normal AND error exit. The temp file
            # name is unique per run, so a run whose lock was reclaimed can never
            # clobber the new owner's in-progress write during cleanup.
            my_owner <- paste(Sys.getpid(), Sys.info()["nodename"], format(Sys.time()))
            my_tmp   <- paste0(reg_path, ".", Sys.getpid(), ".tmp")
            touch_owner <- function() writeLines(my_owner, owner_file)
            touch_owner()
            on.exit({
                # Remove the lock only while it is still OURS: a concurrent run
                # that classified this one as stale may have reclaimed (and
                # replaced) the lock directory under the same path name.
                cur_owner <- if (file.exists(owner_file)) {
                    paste(readLines(owner_file, warn = FALSE), collapse = " ")
                } else ""
                if (nzchar(cur_owner) && cur_owner == my_owner) unlink(lock, recursive = TRUE)
                if (file.exists(my_tmp)) unlink(my_tmp)
            }, add = TRUE)

            # Heartbeat: staleness is judged from owner-file mtime, so touch on
            # wall-clock time — an iteration count is disproportional once the
            # resolution phase or the persist (history copy) gets slow.
            last_touch <- Sys.time()
            heartbeat <- function() {
                now <- Sys.time()
                if (as.numeric(now) - as.numeric(last_touch) > 60) {
                    touch_owner()
                    last_touch <<- now
                }
                invisible(NULL)
            }
            assert_lock_ours <- function() {
                cur_owner <- if (file.exists(owner_file)) {
                    paste(readLines(owner_file, warn = FALSE), collapse = " ")
                } else ""
                if (!nzchar(cur_owner) || cur_owner != my_owner) {
                    stop("catalog registry lock was reclaimed as stale — aborting the write to avoid a corrupt merge")
                }
                invisible(NULL)
            }

            # Read the registry UNDER the lock: the snapshot is fresh, so no
            # stale-write window exists for concurrent runs.
            if (file.exists(reg_path)) {
                rr <- read_registry(reg_path)
                meta <- rr\$meta
                reg_df <- if (!is.null(rr\$data)) normalize_cols(rr\$data) else rr\$data
            } else {
                meta <- list(genome_build = BUILD, reference_annotation = REFANN,
                             catalog_registry_project = PROJ, registry_version = "1")
            }
            if (nn(meta[["genome_build"]]) && meta[["genome_build"]] != BUILD) {
                stop(sprintf("catalog registry built for %s but this run is %s — start a new registry (fail-fast).",
                             meta[["genome_build"]], BUILD))
            }
            if (nn(meta[["reference_annotation"]]) && meta[["reference_annotation"]] != REFANN) {
                stop(sprintf("catalog registry annotated with %s but this run uses %s — start a new registry (fail-fast).",
                             meta[["reference_annotation"]], REFANN))
            }
            # Project provenance guard: the registry carries a semantic owner.
            # Mixing projects into one registry is a provenance risk (silent id
            # reuse across studies), so fail fast unless explicitly allowed.
            meta_proj <- if (is.null(meta[["catalog_registry_project"]])) "" else meta[["catalog_registry_project"]]
            cross_ok <- tolower(as.character("${params.catalog_registry_allow_cross_project}")) %in% c("true", "1", "yes")
            if (nn(meta_proj) && !nn(PROJ) && !cross_ok) {
                stop(sprintf("catalog registry is annotated for project '%s' but this run has no catalog_registry_project — set catalog_registry_project (recommended) or catalog_registry_allow_cross_project=true to extend it",
                             meta_proj))
            }
            if (nn(meta_proj) && nn(PROJ) && meta_proj != PROJ && !cross_ok) {
                stop(sprintf("catalog registry built for project '%s' but this run is '%s' — set catalog_registry_allow_cross_project=true to extend it across projects, or use a separate registry path",
                             meta_proj, PROJ))
            }
            # Backfill the header the first time a project-labelled run touches a
            # legacy registry (rows themselves carry first_seen_project already).
            if (nn(PROJ) && !nn(meta_proj)) meta[["catalog_registry_project"]] <- PROJ

            reg_gene <- if (!is.null(reg_df)) reg_df[reg_df\$entity_type == "gene", ] else data.frame()
            reg_tx   <- if (!is.null(reg_df)) reg_df[reg_df\$entity_type == "transcript", ] else data.frame()
            if (nrow(reg_gene) > 0) global_max <- max(vapply(reg_gene\$catalog_id, gene_num, integer(1)), 0L)
        }

        # ---- Transcript matching against the registry ----
        # Geometry tag from a transcript's exon chain: "S" monoexonic, "M"
        # multiexonic, "U" unknown. The stored signature alone cannot tell a
        # monoexonic exon extent from a 2-exon intron chain (both render as
        # `chr|A-B|strand`), so exact matching keys on signature AND geometry to
        # avoid adopting an unrelated transcript's catalog id.
        geo_of <- function(coords) {
            if (is.na(coords) || !nzchar(coords)) return("U")
            toks <- strsplit(coords, ";", fixed = TRUE)[[1]]
            if (length(toks) >= 2) "M" else "S"
        }
        match_tx <- function(strand, chr, sig, coords) {
            if (is.null(reg_tx) || nrow(reg_tx) == 0) return(list(kind = "none"))
            cand <- reg_tx[reg_tx\$strand == strand & reg_tx\$chr == chr, , drop = FALSE]
            if (nrow(cand) == 0) return(list(kind = "none"))
            sig_ok <- !is.na(cand\$splice_junction_signature) & cand\$splice_junction_signature == sig
            qgeo <- geo_of(coords)
            if (qgeo != "U") {
                cgeo <- vapply(cand\$exon_coords, geo_of, character(1))
                sig_ok <- sig_ok & (cgeo == qgeo | cgeo == "U")
            }
            exact <- cand[sig_ok, , drop = FALSE]
            if (nrow(exact) > 0) {
                u <- unique(exact\$catalog_id)
                if (length(u) == 1) {
                    p <- exact\$parent_catalog_id[1]
                    return(list(kind = "match", catalog_tx_id = u,
                                parent = ifelse(is.na(p) || !nzchar(p), NA_character_, p)))
                }
                return(list(kind = "ambig", candidates = exact))
            }
            r <- vapply(seq_len(nrow(cand)), function(i) exonic_recip(coords, cand\$exon_coords[i]), numeric(1))
            wins <- which(r >= MIN_REC)
            if (length(wins) == 0) return(list(kind = "none"))
            if (length(wins) == 1) {
                w <- cand[wins[1], ]
                return(list(kind = "match", catalog_tx_id = w\$catalog_id,
                            parent = ifelse(is.na(w\$parent_catalog_id) || !nzchar(w\$parent_catalog_id),
                                            NA_character_, w\$parent_catalog_id)))
            }
            list(kind = "ambig", candidates = cand[wins, , drop = FALSE])
        }

        write_review <- function(issue, novel, registry_id, detail) {
            row <- data.frame(issue = issue, novel_gene = novel, registry_id = registry_id,
                              detail = detail, stringsAsFactors = FALSE)
            if (file.exists("registry_review_needed.tsv")) {
                prev <- read_tsv("registry_review_needed.tsv", show_col_types = FALSE)
                row <- suppressWarnings(bind_rows(prev, row))
            }
            write_tsv(row, "registry_review_needed.tsv")
        }

        # ---- Resolve catalog IDs gene by gene (corrections A & B) ----
        gene_order <- tx_dat %>%
            group_by(gene_id) %>%
            summarise(chr = dplyr::first(chr), s = min(span_start), .groups = "drop") %>%
            arrange(chr, s)
        gene_ids <- gene_order\$gene_id

        next_gene_num <- global_max + 1L
        tx_map <- list()          # novel transcript_id -> catalog_tx_id

        catalog_ids_out <- data.frame(
            gene_id = character(), mstrg_gene_id = character(), catalog_id = character(),
            locus = character(), chr = character(), start = integer(), end = integer(),
            strand = character(), class_codes = character(), stringsAsFactors = FALSE)
        catalog_tx_out <- data.frame(
            catalog_tx_id = character(), mstrg_tx_id = character(),
            parent_catalog_id = character(), firma = character(), stringsAsFactors = FALSE)
        new_gene_rows <- data.frame()
        new_tx_rows   <- data.frame()

        # Heartbeat: keeps the lock owner fresh while resolving (time-based).
        for (gid in gene_ids) {
            heartbeat()
            idx <- which(tx_dat\$gene_id == gid)
            gdat <- tx_dat[idx, ]
            gchr <- gdat\$chr[1]; gstrand <- gdat\$strand[1]
            gstart <- min(gdat\$span_start); gend <- max(gdat\$span_end)
            class_codes <- paste(sort(unique(gdat\$class_code[!is.na(gdat\$class_code)])), collapse = ",")
            if (!nzchar(class_codes)) class_codes <- NA_character_

            sig_rows <- unique(gdat[, c("signature", "exon_coords")])
            matched_tx <- list()
            used_registry_tx <- character()
            inherited_parents <- character()

            for (k in seq_len(nrow(sig_rows))) {
                sg  <- sig_rows\$signature[k]; cds <- sig_rows\$exon_coords[k]
                gkey <- paste0(geo_of(cds), "|", sg)
                res <- match_tx(gstrand, gchr, sg, cds)
                if (res\$kind == "ambig") {
                    for (ci in seq_len(nrow(res\$candidates))) {
                        write_review("transcript ambiguity",
                                     paste0(gid, " :: ", sg),
                                     res\$candidates\$catalog_id[ci],
                                     paste("exonic/structural overlap >= min_reciprocal; different entities"))
                    }
                    stop("AMBIGUOUS transcript match for ", gid, " (signature ", sg,
                         "). Inspect registry_review_needed.tsv and reconcile the registry manually.")
                }
                if (res\$kind == "match") {
                    matched_tx[[gkey]] <- res\$catalog_tx_id
                    used_registry_tx <- c(used_registry_tx, res\$catalog_tx_id)
                    inherited_parents <- c(inherited_parents, res\$parent)
                }
            }
            if (length(used_registry_tx) != length(unique(used_registry_tx))) {
                write_review("signature collision", gid, paste(used_registry_tx, collapse = ","),
                             "distinct isoforms of the same gene matched the same registry transcript")
                stop("Signature collision in ", gid, ": two isoforms map to the same catalog_tx — inspect registry_review_needed.tsv.")
            }
            nonNA_parents <- unique(inherited_parents[!is.na(inherited_parents)])

            if (length(nonNA_parents) >= 2) {
                write_review("parent ambiguity", gid, paste(nonNA_parents, collapse = ","),
                             "transcripts of the same gene inherited distinct parent_catalog_id")
                stop("AMBIGUOUS parent catalog for ", gid, ": inherited parents = ",
                     paste(nonNA_parents, collapse = ","), " — inspect registry_review_needed.tsv.")
            }

            if (length(nonNA_parents) == 1) {
                gene_catalog_id <- nonNA_parents[1]
                gene_from <- "inherited"
            } else {
                # No transcript matched: gene-level span fallback, else a fresh id
                gene_catalog_id <- NA_character_; gene_from <- "none"
                if (!reg_off && nrow(reg_gene) > 0) {
                    gcand <- reg_gene[reg_gene\$strand == gstrand & reg_gene\$chr == gchr, , drop = FALSE]
                    if (nrow(gcand) > 0) {
                        r <- vapply(seq_len(nrow(gcand)), function(i)
                            span_recip(gstart, gend, gcand\$start[i], gcand\$end[i]), numeric(1))
                        wins <- which(r >= MIN_REC)
                        if (length(wins) >= 2) {
                            for (wi in wins) write_review("gene span ambiguity", gid, gcand\$catalog_id[wi],
                                                          "gene span overlaps >= min_reciprocal with multiple registry genes")
                            stop("AMBIGUOUS gene span match for ", gid, " — inspect registry_review_needed.tsv.")
                        }
                        if (length(wins) == 1) { gene_catalog_id <- gcand\$catalog_id[wins[1]]; gene_from <- "span" }
                    }
                }
                if (is.na(gene_catalog_id)) {
                    gene_catalog_id <- sprintf("%s%06d", PREFIX, next_gene_num)
                    next_gene_num <- next_gene_num + 1L
                    gene_from <- "new"
                }
            }

            # Next free transcript ordinal under this gene (reused when the gene is inherited)
            existing_suffixes <- integer()
            if (!reg_off && nrow(reg_tx) > 0) {
                keep_tx <- which(reg_tx\$parent_catalog_id == gene_catalog_id)
                if (length(keep_tx) > 0) {
                    existing_suffixes <- as.integer(vapply(reg_tx\$catalog_id[keep_tx], tx_suffix, integer(1)))
                    existing_suffixes <- existing_suffixes[!is.na(existing_suffixes)]
                }
            }
            next_tx_num <- if (length(existing_suffixes) > 0) max(existing_suffixes) + 1L else 1L

            # Per-gene assignment keyed by signature (new ids get the next free ordinal)
            gene_assign <- list()
            for (sg in names(matched_tx)) gene_assign[[sg]] <- matched_tx[[sg]]

            g_order <- gdat\$transcript_id[order(gdat\$span_start, gdat\$span_end)]
            for (tid in g_order) {
                r <- which(gdat\$transcript_id == tid)
                sg  <- gdat\$signature[r]; cds <- gdat\$exon_coords[r]
                cc  <- gdat\$class_code[r]; st <- gdat\$span_start[r]; en <- gdat\$span_end[r]
                gkey <- paste0(geo_of(cds), "|", sg)
                if (!(gkey %in% names(gene_assign))) {
                    cx <- sprintf("%s.%d", gene_catalog_id, next_tx_num)
                    next_tx_num <- next_tx_num + 1L
                    gene_assign[[gkey]] <- cx
                    new_tx_rows <- rbind(new_tx_rows, data.frame(
                        entity_type = "transcript", catalog_id = cx,
                        parent_catalog_id = gene_catalog_id, mstrg_id = tid,
                        chr = gchr, start = st, end = en, strand = gstrand, class_code = cc,
                        splice_junction_signature = sg, exon_coords = cds,
                        first_seen_project = PROJ, first_seen_run = RUN_TOKEN,
                        first_seen_date = format(Sys.time()), stringsAsFactors = FALSE))
                }
                cx_id <- gene_assign[[gkey]]
                tx_map[[tid]] <- c(catalogs = cx_id, gene = gene_catalog_id, old_gene = gid)
                catalog_tx_out <- rbind(catalog_tx_out, data.frame(
                    catalog_tx_id = cx_id, mstrg_tx_id = tid,
                    parent_catalog_id = gene_catalog_id, firma = sg,
                    stringsAsFactors = FALSE))
            }

            if (gene_from == "new") {
                new_gene_rows <- rbind(new_gene_rows, data.frame(
                    entity_type = "gene", catalog_id = gene_catalog_id,
                    parent_catalog_id = NA_character_, mstrg_id = gid,
                    chr = gchr, start = gstart, end = gend, strand = gstrand,
                    class_code = class_codes, splice_junction_signature = NA_character_,
                    exon_coords = NA_character_, first_seen_project = PROJ,
                    first_seen_run = RUN_TOKEN, first_seen_date = format(Sys.time()),
                    stringsAsFactors = FALSE))
            }

            locus <- sprintf("%s:%d-%d(%s)", gchr, gstart, gend, gstrand)
            catalog_ids_out <- rbind(catalog_ids_out, data.frame(
                gene_id = gid, mstrg_gene_id = gid, catalog_id = gene_catalog_id,
                locus = locus, chr = gchr, start = gstart, end = gend,
                strand = gstrand, class_codes = class_codes, stringsAsFactors = FALSE))
        }

        # ---- Assertions (plan T2.6) ----
        if (nrow(catalog_ids_out) != length(unique(catalog_ids_out\$gene_id)))
            stop("Each novel gene must map to exactly one catalog_id — duplication detected")
        if (any(duplicated(catalog_ids_out\$catalog_id)))
            stop("Duplicate gene catalog_id assigned")
        if (nrow(catalog_tx_out) != length(unique(catalog_tx_out\$catalog_tx_id)))
            stop("Each novel transcript must map to exactly one catalog_tx_id — duplication detected")
        if (length(tx_map) != nrow(tx_dat))
            stop("Every parsed frozen transcript must receive a catalog_tx_id")

        # ---- Emit catalog tables ----
        cat("Assigned - novel genes:", nrow(catalog_ids_out),
            "- novel transcripts:", nrow(catalog_tx_out), "\\n")
        write_tsv(catalog_ids_out, "catalog_ids.tsv")
        write_tsv(catalog_tx_out, "catalog_tx.tsv")

        # ---- Regenerate frozen GTF with catalog IDs (export artifact only) ----
        out_lines <- gtf_lines
        for (tid in names(tx_map)) {
            cx <- tx_map[[tid]]["catalogs"]; ng <- tx_map[[tid]]["gene"]; og <- tx_map[[tid]]["old_gene"]
            pat_tx <- paste0('transcript_id "', tid, '"')
            if (!any(grepl(pat_tx, out_lines, fixed = TRUE))) next
            pat_ge <- paste0('gene_id "', og, '"')
            rep_ge <- paste0('gene_id "', ng, '"')
            rep_tx <- paste0('transcript_id "', cx, '"')
            hit <- grepl(pat_tx, out_lines, fixed = TRUE)
            out_lines[hit] <- gsub(pat_ge, rep_ge, out_lines[hit], fixed = TRUE)
            out_lines[hit] <- gsub(pat_tx, rep_tx, out_lines[hit], fixed = TRUE)
            out_lines[hit] <- paste0(out_lines[hit],
                sprintf(' mstrg_gene_id "%s"; mstrg_transcript_id "%s"; catalog_id "%s";', og, tid, ng))
        }
        writeLines(out_lines, "frozen_lncrna.catalog.gtf")
        cat("Wrote frozen_lncrna.catalog.gtf with", nrow(catalog_tx_out), "cataloged transcripts\\n")

        # ---- Persist registry (atomic, with history backup) ----
        if (!reg_off) {
            all_new <- rbind(new_gene_rows, new_tx_rows)
            if (nrow(all_new) > 0) all_new <- normalize_cols(all_new)
            # Skip the rewrite + history backup when nothing changed and the
            # on-disk file is already canonical (one-time heal of legacy files
            # with trailing whitespace still triggers a rewrite).
            needs_rewrite <- nrow(all_new) > 0
            if (!needs_rewrite && file.exists(reg_path)) {
                cur <- tryCatch(readLines(reg_path, warn = FALSE), error = function(e) character(0))
                if (any(endsWith(cur, " "))) needs_rewrite <- TRUE
            }
            if (!needs_rewrite) {
                cat("Registry unchanged (no new entities, schema canonical):", reg_path, "\\n")
            } else {
                reg_out <- if (!is.null(reg_df) && nrow(reg_df) > 0) {
                    if (nrow(all_new) > 0) rbind(reg_df, all_new) else reg_df
                } else {
                    all_new
                }
                if (any(duplicated(reg_out\$catalog_id[reg_out\$entity_type == "gene"])))
                    stop("Registry integrity: duplicate gene catalog_id would be written")
                if (any(duplicated(reg_out\$catalog_id[reg_out\$entity_type == "transcript"])))
                    stop("Registry integrity: duplicate transcript catalog_id would be written")
                # The lock could have been reclaimed mid-persist by a stale
                # decision; never write into a registry we no longer own.
                heartbeat()
                assert_lock_ours()
                tmp_path <- my_tmp
                if (file.exists(tmp_path)) unlink(tmp_path)
                write_registry(tmp_path, meta, reg_out)
                if (file.exists(reg_path)) {
                    bk <- file.path(rdir, "history",
                                    sprintf("lncrna_catalog_registry.%s.%s.tsv",
                                            format(Sys.time(), "%Y%m%dT%H%M%S"), Sys.getpid()))
                    if (!file.copy(reg_path, bk)) stop("Could not back up previous registry to ", bk)
                }
                assert_lock_ours()
                if (!file.rename(tmp_path, reg_path)) {
                    stop("Could not atomically move registry into place: ", tmp_path, " -> ", reg_path)
                }
                cat("Registry updated:", reg_path, "-", nrow(reg_out), "rows\\n")
            }
        }
    }

    run_registry_write()
    """
}