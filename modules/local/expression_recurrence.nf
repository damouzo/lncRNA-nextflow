// Expression recurrence filter (condition-blind)
// Reads featureCounts matrices from all samples, computes per-transcript
// expression prevalence, and keeps transcripts expressed in >= N replicates

process EXPRESSION_RECURRENCE {

    input:
    path catalog_table        // from catalog_overlap
    path fc_tsv               // collected featureCounts TSVs from all samples
    path merged_gtf           // merged StringTie2 GTF (for pass-through)

    output:
    path "recurrence_filtered.tsv", emit: recurrence_table
    path "recurrence_filtered.gtf", emit: recurrence_gtf
    path "recurrence_stats.txt", emit: stats

    script:
    """
    #!/usr/bin/env Rscript
    library(dplyr)
    library(readr)
    library(tibble)

    MIN_REPS <- ${params.expression_recurrence_min_reps}
    CPM_THRESHOLD <- 1.0

    extract_attr <- function(attrs, key) {
        pat <- paste0('(?:^|; )', key, ' "([^"]+)"')
        m <- regexpr(pat, attrs, perl = TRUE)
        ans <- rep(NA_character_, length(attrs))
        ok <- m > 0L
        if (any(ok)) ans[ok] <- sub(paste0('.*', pat, '.*\$'), '\\\\1', attrs[ok], perl = TRUE)
        ans
    }

    # Parse class_code from merged GTF by transcript_id
    gtf_lines <- readLines("${merged_gtf}")
    gtf_parts <- strsplit(gtf_lines, "\\t", fixed = TRUE)
    gtf_ok <- lengths(gtf_parts) >= 9L
    gtf_attrs <- vapply(gtf_parts[gtf_ok], function(p) p[[9L]], character(1))
    gtf_has_tx <- grepl("transcript_id", gtf_attrs, fixed = TRUE)
    tx_cc <- data.frame(
        transcript_id = extract_attr(gtf_attrs[gtf_has_tx], "transcript_id"),
        class_code    = extract_attr(gtf_attrs[gtf_has_tx], "class_code"),
        stringsAsFactors = FALSE
    ) %>% filter(!is.na(transcript_id), !is.na(class_code)) %>%
        distinct(transcript_id, .keep_all = TRUE)

    # Read the catalog table
    catalog <- read_tsv("${catalog_table}", show_col_types = FALSE)

    # Find and read all featureCounts output files
    # Nextflow stages them in the work directory
    fc_files <- list.files(pattern = "featureCounts", full.names = TRUE)

    if (length(fc_files) == 0) {
        cat("WARNING: no featureCounts files found, skipping recurrence filter\\n")
        catalog\$n_reps_expressed <- NA_integer_
        catalog\$passes_recurrence <- TRUE
        write_tsv(catalog, "recurrence_filtered.tsv")
        file.copy("${merged_gtf}", "recurrence_filtered.gtf")
        quit(save = "no", status = 0)
    }

    cat("Found", length(fc_files), "featureCounts files\\n")

    # Parse each featureCounts file and build count matrix
    count_list <- list()
    sample_names <- character()

    for (f in fc_files) {
        # featureCounts output: skip comment lines (starting with #), first line is header
        # Columns: Geneid, Chr, Start, End, Strand, Length, <bam_path>
        lines <- readLines(f)
        data_lines <- lines[!grepl("^#", lines) & !grepl("^\$", lines)]

        if (length(data_lines) < 2) next  # only header, no data

        header <- strsplit(data_lines[1], "\\t")[[1]]
        # The last column is the count column (sample name)
        count_col <- length(header)

        df <- read_tsv(f, comment = "#", show_col_types = FALSE, col_names = FALSE, skip = 2)

        if (nrow(df) == 0) next

        # featureCounts adds summary lines at the bottom; exclude them
        # Summary rows have Geneid starting with non-gene-like names
        df <- df[!grepl("^(Assigned|Unassigned|__)", df[[1]]), , drop = FALSE]

        tx_ids <- df[[1]]
        counts <- as.numeric(df[[count_col]])
        names(counts) <- tx_ids

        count_list[[length(count_list) + 1]] <- counts

        # Extract sample name from filename (before .featureCounts.tsv)
        sn <- sub(".featureCounts.tsv", "", basename(f), fixed = TRUE)
        sample_names <- c(sample_names, sn)
    }

    if (length(count_list) == 0) {
        cat("WARNING: could not parse any featureCounts files\\n")
        catalog\$n_reps_expressed <- NA_integer_
        catalog\$passes_recurrence <- TRUE
        write_tsv(catalog, "recurrence_filtered.tsv")
        file.copy("${merged_gtf}", "recurrence_filtered.gtf")
        quit(save = "no", status = 0)
    }

    # Build count matrix: rows = transcripts, columns = samples
    all_tx <- unique(unlist(lapply(count_list, names)))
    count_matrix <- matrix(0, nrow = length(all_tx), ncol = length(sample_names))
    rownames(count_matrix) <- all_tx
    colnames(count_matrix) <- sample_names

    for (i in seq_along(count_list)) {
        tx <- names(count_list[[i]])
        count_matrix[tx, i] <- count_list[[i]][tx]
    }

    cat(sprintf("Count matrix: %d transcripts x %d samples\\n",
                nrow(count_matrix), ncol(count_matrix)))

    # Compute CPM and count samples with expression above threshold
    lib_sizes <- colSums(count_matrix)
    cpm <- sweep(count_matrix, 2, lib_sizes / 1e6, "/")
    expressed <- cpm > CPM_THRESHOLD
    n_expressed <- rowSums(expressed)

    # Build recurrence table
    rec_df <- data.frame(
        transcript_id = rownames(count_matrix),
        n_reps_expressed = n_expressed,
        total_samples = ncol(count_matrix),
        max_cpm = apply(cpm, 1, max),
        mean_cpm = rowMeans(cpm),
        stringsAsFactors = FALSE
    )

    # Match to catalog table by transcript ID
    # The catalog may use different ID column names
    id_col <- "ID"
    if (!(id_col %in% colnames(catalog))) {
        # Try common alternatives
        for (cname in colnames(catalog)) {
            if (grepl("transcript_id|ID|gene_id", cname, ignore.case = TRUE)) {
                id_col <- cname
                break
            }
        }
    }

    catalog\$n_reps_expressed <- rec_df\$n_reps_expressed[
        match(catalog[[id_col]], rec_df\$transcript_id)]
    catalog\$n_reps_expressed[is.na(catalog\$n_reps_expressed)] <- 0
    catalog\$passes_recurrence <- catalog\$n_reps_expressed >= MIN_REPS

    catalog\$class_code <- tx_cc\$class_code[match(catalog[[id_col]], tx_cc\$transcript_id)]

    n_pass <- sum(catalog\$passes_recurrence, na.rm = TRUE)
    n_total <- nrow(catalog)
    cat(sprintf("Recurrence filter: %d / %d transcripts pass (>= %d reps @ CPM > %.1f)\\n",
                n_pass, n_total, MIN_REPS, CPM_THRESHOLD))

    write_tsv(catalog, "recurrence_filtered.tsv")
    file.copy("${merged_gtf}", "recurrence_filtered.gtf")

    # Summary stats
    sink("recurrence_stats.txt")
    cat(sprintf("total_transcripts\\t%d\\n", n_total))
    cat(sprintf("passed_recurrence\\t%d\\n", n_pass))
    cat(sprintf("min_reps_threshold\\t%d\\n", MIN_REPS))
    cat(sprintf("cpm_threshold\\t%.1f\\n", CPM_THRESHOLD))
    cat("n_reps_distribution:\\n")
    print(table(catalog\$n_reps_expressed))
    sink()

    cat("DONE\\n")
    """
}
