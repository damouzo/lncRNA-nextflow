#!/usr/bin/env Rscript
# Expression recurrence filter (condition-blind). Keeps transcripts expressed in
# >= min-reps samples above a CPM threshold, where CPM is normalised by the total
# fragments of each sample (every category of the featureCounts .summary), never
# by the candidate-only column sums.
#
# usage: expression_recurrence.R --catalog <tsv> --merged-gtf <gtf> \
#          --fc-dir <dir> --min-reps 3 --cpm-threshold 1.0 --outdir <dir>

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tibble)
})

get_arg <- function(flag, default = NULL) {
    args <- commandArgs(trailingOnly = TRUE)
    i <- match(flag, args)
    if (is.na(i)) return(default)
    args[i + 1]
}

catalog_table <- get_arg("--catalog")
merged_gtf    <- get_arg("--merged-gtf")
fc_dir        <- get_arg("--fc-dir", ".")
MIN_REPS      <- as.integer(get_arg("--min-reps", "3"))
CPM_THRESHOLD <- as.numeric(get_arg("--cpm-threshold", "1.0"))
outdir        <- get_arg("--outdir", ".")
if (is.null(catalog_table) || is.null(merged_gtf)) {
    stop("--catalog and --merged-gtf are required")
}
out <- function(f) file.path(outdir, f)

extract_attr <- function(attrs, key) {
    pat <- paste0('(?:^|; )', key, ' "([^"]+)"')
    m <- regexpr(pat, attrs, perl = TRUE)
    ans <- rep(NA_character_, length(attrs))
    ok <- m > 0L
    if (any(ok)) ans[ok] <- sub(paste0('.*', pat, '.*$'), '\\1', attrs[ok], perl = TRUE)
    ans
}

# Parse class_code from merged GTF by transcript_id
gtf_lines <- readLines(merged_gtf)
gtf_parts <- strsplit(gtf_lines, "\t", fixed = TRUE)
gtf_ok <- lengths(gtf_parts) >= 9L
gtf_attrs <- vapply(gtf_parts[gtf_ok], function(p) p[[9L]], character(1))
gtf_has_tx <- grepl("transcript_id", gtf_attrs, fixed = TRUE)
tx_cc <- data.frame(
    transcript_id = extract_attr(gtf_attrs[gtf_has_tx], "transcript_id"),
    class_code    = extract_attr(gtf_attrs[gtf_has_tx], "class_code"),
    stringsAsFactors = FALSE
) %>% filter(!is.na(transcript_id), !is.na(class_code)) %>%
    distinct(transcript_id, .keep_all = TRUE)

catalog <- read_tsv(catalog_table, show_col_types = FALSE)

id_col <- "ID"
if (!(id_col %in% colnames(catalog))) {
    for (cname in colnames(catalog)) {
        if (grepl("transcript_id|ID|gene_id", cname, ignore.case = TRUE)) {
            id_col <- cname
            break
        }
    }
}
# Set before any early exit: ANNOTATION_FREEZE gates on class_code.
catalog$class_code <- tx_cc$class_code[match(catalog[[id_col]], tx_cc$transcript_id)]

# Count files only: exclude the .summary files (they share the "featureCounts" stem)
fc_files <- list.files(path = fc_dir, pattern = "featureCounts", full.names = TRUE)
fc_files <- fc_files[!grepl(".summary", fc_files, fixed = TRUE)]

if (length(fc_files) == 0) {
    cat("WARNING: no featureCounts files found, skipping recurrence filter\n")
    catalog$n_reps_expressed <- NA_integer_
    catalog$passes_recurrence <- TRUE
    write_tsv(catalog, out("recurrence_filtered.tsv"))
    file.copy(merged_gtf, out("recurrence_filtered.gtf"))
    quit(save = "no", status = 0)
}

cat("Found", length(fc_files), "featureCounts files\n")

count_list <- list()
sample_names <- character()

for (f in fc_files) {
    data_lines <- readLines(f)
    data_lines <- data_lines[!grepl("^#", data_lines) & !grepl("^$", data_lines)]
    if (length(data_lines) < 2) next

    header <- strsplit(data_lines[1], "\t")[[1]]
    count_col <- length(header)

    df <- read_tsv(f, comment = "#", show_col_types = FALSE, col_names = FALSE, skip = 2)
    if (nrow(df) == 0) next

    df <- df[!grepl("^(Assigned|Unassigned|__)", df[[1]]), , drop = FALSE]
    if (nrow(df) == 0) next

    counts <- as.numeric(df[[count_col]])
    names(counts) <- df[[1]]
    count_list[[length(count_list) + 1]] <- counts
    sample_names <- c(sample_names, sub(".featureCounts.tsv", "", basename(f), fixed = TRUE))
}

if (length(count_list) == 0) {
    cat("WARNING: could not parse any featureCounts files\n")
    catalog$n_reps_expressed <- NA_integer_
    catalog$passes_recurrence <- TRUE
    write_tsv(catalog, out("recurrence_filtered.tsv"))
    file.copy(merged_gtf, out("recurrence_filtered.gtf"))
    quit(save = "no", status = 0)
}

all_tx <- unique(unlist(lapply(count_list, names)))
count_matrix <- matrix(0, nrow = length(all_tx), ncol = length(sample_names))
rownames(count_matrix) <- all_tx
colnames(count_matrix) <- sample_names
for (i in seq_along(count_list)) {
    tx <- names(count_list[[i]])
    count_matrix[tx, i] <- count_list[[i]][tx]
}

cat(sprintf("Count matrix: %d transcripts x %d samples\n",
            nrow(count_matrix), ncol(count_matrix)))

# H2 / 3.1: CPM denominator = total fragments of the sample = every category of the
# featureCounts .summary (Assigned + all Unassigned_*). The candidate-only column
# sums ignore reads outside the candidate set and would inflate CPM.
all_files <- list.files()
summary_files <- all_files[grepl(".featureCounts.tsv.summary", all_files, fixed = TRUE)]
if (length(summary_files) > 0) {
    tot <- numeric(0)
    for (i in seq_along(summary_files)) {
        f <- summary_files[[i]]
        sn <- sub(".featureCounts.tsv.summary", "", basename(f), fixed = TRUE)
        tbl <- read.delim(f, header = FALSE, stringsAsFactors = FALSE,
                          comment.char = "", quote = "")
        vals <- suppressWarnings(as.numeric(tbl[[2]][-1]))
        tot[sn] <- sum(vals, na.rm = TRUE)
    }
    lib_sizes <- tot[colnames(count_matrix)]
    if (any(is.na(lib_sizes)) || any(lib_sizes <= 0)) {
        cat("WARNING: featureCounts .summary incomplete; falling back to candidate column sums\n")
        lib_sizes <- colSums(count_matrix)
    } else {
        cat(sprintf("CPM denominator: per-sample total fragments (min %g, max %g)\n",
                    min(lib_sizes), max(lib_sizes)))
    }
} else {
    cat("WARNING: no featureCounts .summary files; falling back to candidate column sums\n")
    lib_sizes <- colSums(count_matrix)
}
cpm <- sweep(count_matrix, 2, lib_sizes / 1e6, "/")
expressed <- cpm > CPM_THRESHOLD
n_expressed <- rowSums(expressed)

rec_df <- data.frame(
    transcript_id = rownames(count_matrix),
    n_reps_expressed = n_expressed,
    total_samples = ncol(count_matrix),
    max_cpm = apply(cpm, 1, max),
    mean_cpm = rowMeans(cpm),
    stringsAsFactors = FALSE
)

catalog$n_reps_expressed <- rec_df$n_reps_expressed[
    match(catalog[[id_col]], rec_df$transcript_id)]
catalog$n_reps_expressed[is.na(catalog$n_reps_expressed)] <- 0
catalog$passes_recurrence <- catalog$n_reps_expressed >= MIN_REPS

n_pass <- sum(catalog$passes_recurrence, na.rm = TRUE)
n_total <- nrow(catalog)
cat(sprintf("Recurrence filter: %d / %d transcripts pass (>= %d reps @ CPM > %.1f)\n",
            n_pass, n_total, MIN_REPS, CPM_THRESHOLD))

write_tsv(catalog, out("recurrence_filtered.tsv"))
file.copy(merged_gtf, out("recurrence_filtered.gtf"))

sink(out("recurrence_stats.txt"))
cat(sprintf("total_transcripts\t%d\n", n_total))
cat(sprintf("passed_recurrence\t%d\n", n_pass))
cat(sprintf("min_reps_threshold\t%d\n", MIN_REPS))
cat(sprintf("cpm_threshold\t%.1f\n", CPM_THRESHOLD))
cat("cpm_denominator\tper_sample_total_fragments\n")
cat("n_reps_distribution:\n")
print(table(catalog$n_reps_expressed))
sink()

cat("DONE\n")
