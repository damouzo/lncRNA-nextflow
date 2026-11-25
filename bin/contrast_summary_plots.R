#!/usr/bin/env Rscript
# Cross-contrast comparison figures + tables, split by lncRNA subset
# (combined = annotated + novel / annotated / novel), mirroring the volcano tabs.
# Writes everything under <outdir>/plots/<subset>/ and <outdir>/tables/:
#   plots/<subset>/sig_counts_by_contrast.{png,pdf}   significant up/down per contrast
#   plots/<subset>/venn_lncrna.{png,pdf}              Venn of significant lncRNAs (<= 4 contrasts)
#   plots/<subset>/upset_lncrna.{png,pdf}             membership overlap (any contrast count)
#   tables/membership_<subset>.tsv                    gene x contrast 0/1 membership
#   tables/lncrna_across_contrasts.tsv                wide copy + significant_contrasts column
#
# Non-blocking: every panel is guarded so a missing package or an empty input
# never fails the task (the module also sets errorStrategy 'ignore'). Each panel
# is written as PNG (slides) and PDF (print/vector).
#
#   contrast_summary_plots.R --across lncrna_across_contrasts.tsv --outdir .

suppressPackageStartupMessages({
  library(readr)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}
across_path <- getarg("--across")
outdir      <- getarg("--outdir", ".")
fdr         <- as.numeric(getarg("--fdr", "0.05"))
plots_dir   <- file.path(outdir, "plots")
tables_dir  <- file.path(outdir, "tables")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# ggvenn / ComplexUpset are extras: guard their absence so it only skips the
# Venn/UpSet panels, never the whole task.
have <- function(pkg) requireNamespace(pkg, quietly = TRUE)

SUBSETS <- list(combined = c("novel", "annotated_lncrna"),
                annotated = "annotated_lncrna",
                novel = "novel")
# ggvenn reads 2-4 sets; beyond that a Venn diagram is not interpretable.
VENN_MAX_SETS <- 4L

save_panel <- function(subset, name, plot, width, height) {
  d <- file.path(plots_dir, subset)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  tryCatch(ggplot2::ggsave(file.path(d, paste0(name, ".png")), plot,
                           width = width, height = height, dpi = 150),
           error = function(e) cat(sprintf("%s/%s PNG skipped: %s\n", subset, name, conditionMessage(e))))
  tryCatch(ggplot2::ggsave(file.path(d, paste0(name, ".pdf")), plot,
                           width = width, height = height),
           error = function(e) cat(sprintf("%s/%s PDF skipped: %s\n", subset, name, conditionMessage(e))))
}

text_panel <- function(msg) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = msg, size = 4.5, hjust = 0.5, vjust = 0.5) +
    theme_void()
}

ac <- tryCatch(read_tsv(across_path, show_col_types = FALSE), error = function(e) NULL)
if (is.null(ac) || nrow(ac) == 0) {
  cat("contrast_summary_plots: across table empty or unreadable; nothing to do\n")
  quit(save = "no", status = 0)
}

padj_cols <- grep("\\.padj$", colnames(ac), value = TRUE)
ct_ids    <- sub("\\.padj$", "", padj_cols)
if (length(ct_ids) < 1) {
  cat("contrast_summary_plots: no per-contrast .padj columns\n")
  quit(save = "no", status = 0)
}

origin <- if ("origin" %in% colnames(ac)) as.character(ac$origin) else rep(NA_character_, nrow(ac))

sig_of <- function(cid) {
  col <- paste0(cid, ".padj")
  v <- if (col %in% colnames(ac)) suppressWarnings(as.numeric(ac[[col]])) else rep(NA_real_, nrow(ac))
  !is.na(v) & v < fdr
}
lfc_of <- function(cid) {
  col <- paste0(cid, ".log2FoldChange")
  if (!col %in% colnames(ac)) return(rep(NA_real_, nrow(ac)))
  suppressWarnings(as.numeric(ac[[col]]))
}

# ---- per-subset: significant counts, membership table, Venn, UpSet ----
empty_membership <- function(ct_ids) {
  tab <- data.frame(gene_id = character(), gene_id_stable = character(),
                    gene_name = character(), origin = character(),
                    stringsAsFactors = FALSE)
  for (cid in ct_ids) tab[[cid]] <- integer()
  tab$n_contrasts_significant <- integer()
  tab
}

for (sub in names(SUBSETS)) {
  in_sub <- !is.na(origin) & origin %in% SUBSETS[[sub]]
  if (!any(in_sub)) {
    cat(sprintf("subset %s: no lncRNAs; writing empty membership, skipping plots\n", sub))
    tryCatch(write_tsv(empty_membership(ct_ids),
                       file.path(tables_dir, sprintf("membership_%s.tsv", sub))),
             error = function(e) cat(sprintf("membership_%s.tsv skipped: %s\n", sub, conditionMessage(e))))
    next
  }

  memb <- data.frame(gene_id = ac$gene_id[in_sub], stringsAsFactors = FALSE)
  lfc  <- data.frame(gene_id = ac$gene_id[in_sub], stringsAsFactors = FALSE)
  for (cid in ct_ids) {
    memb[[cid]] <- as.integer(sig_of(cid)[in_sub])
    lfc[[cid]]  <- lfc_of(cid)[in_sub]
  }
  member_mat <- as.matrix(memb[, ct_ids, drop = FALSE])

  # --- membership table (query-friendly: one row per gene, 0/1 per contrast) ---
  tryCatch({
    tab <- data.frame(
      gene_id = ac$gene_id[in_sub],
      gene_id_stable = if ("gene_id_stable" %in% colnames(ac)) ac$gene_id_stable[in_sub] else ac$gene_id[in_sub],
      gene_name = if ("gene_name" %in% colnames(ac)) ac$gene_name[in_sub] else NA_character_,
      origin = origin[in_sub],
      stringsAsFactors = FALSE)
    for (cid in ct_ids) tab[[cid]] <- memb[[cid]]
    tab$n_contrasts_significant <- rowSums(member_mat)
    write_tsv(tab, file.path(tables_dir, sprintf("membership_%s.tsv", sub)))
  }, error = function(e) cat(sprintf("membership_%s.tsv skipped: %s\n", sub, conditionMessage(e))))

  # --- significant counts (up/down) per contrast ---
  tryCatch({
    long <- do.call(rbind, lapply(ct_ids, function(cid) {
      s <- memb[[cid]] == 1L
      f <- lfc[[cid]]
      data.frame(contrast = cid,
                 up = sum(s & !is.na(f) & f > 0),
                 down = sum(s & !is.na(f) & f < 0),
                 stringsAsFactors = FALSE)
    }))
    llong <- long %>%
      pivot_longer(c(up, down), names_to = "dir", values_to = "n") %>%
      mutate(contrast = sub("__", " / ", contrast))
    if (sum(llong$n) > 0) {
      p <- ggplot(llong, aes(reorder(contrast, n, FUN = sum), n, fill = dir)) +
        geom_col(position = "dodge") +
        scale_fill_manual(values = c(up = "#b71c1c", down = "#1565c0")) +
        coord_flip() +
        labs(title = sprintf("Significant lncRNAs per contrast — %s", sub),
             x = NULL, y = "N significant (padj < FDR)", fill = NULL) +
        theme_minimal()
      save_panel(sub, "sig_counts_by_contrast", p, 9, max(4, 0.4 * length(ct_ids) + 2))
    } else {
      save_panel(sub, "sig_counts_by_contrast",
                 text_panel("No significant lncRNAs in any contrast for this subset."), 9, 4)
    }
  }, error = function(e) cat(sprintf("%s sig_counts skipped: %s\n", sub, conditionMessage(e))))

  # --- Venn (skipped with a text panel above the interpretable set limit) ---
  tryCatch({
    if (length(ct_ids) > VENN_MAX_SETS) {
      save_panel(sub, "venn_lncrna",
                 text_panel(sprintf(
                   "Too many contrasts (%d): exceeds the interpretable limit of a Venn diagram (max %d sets).",
                   length(ct_ids), VENN_MAX_SETS)), 8, 6)
    } else if (!have("ggvenn")) {
      cat("ggvenn not available; venn_lncrna skipped\n")
    } else {
      sets <- lapply(ct_ids, function(cid) memb$gene_id[memb[[cid]] == 1L])
      names(sets) <- sub("__", " / ", ct_ids)
      sets <- sets[vapply(sets, length, integer(1)) > 0]
      if (length(sets) >= 2) {
        save_panel(sub, "venn_lncrna", ggvenn::ggvenn(sets, show_percentage = FALSE), 8, 6)
      } else {
        save_panel(sub, "venn_lncrna",
                   text_panel("Fewer than two contrasts have significant lncRNAs."), 8, 6)
      }
    }
  }, error = function(e) cat(sprintf("%s venn skipped: %s\n", sub, conditionMessage(e))))

  # --- UpSet (scales to any number of contrasts) ---
  tryCatch({
    bin <- as.data.frame(lapply(ct_ids, function(cid) memb[[cid]] == 1L))
    colnames(bin) <- ct_ids
    bin <- bin[, colSums(bin) > 0, drop = FALSE]
    if (!have("ComplexUpset")) {
      cat("ComplexUpset not available; upset_lncrna skipped\n")
    } else if (ncol(bin) >= 2) {
      bin$gene_id <- memb$gene_id  # extra column, ignored by intersect
      p <- ComplexUpset::upset(bin, intersect = colnames(bin)[colnames(bin) != "gene_id"],
                               sort_sets = FALSE)
      save_panel(sub, "upset_lncrna", p, 12, 7)
    } else {
      save_panel(sub, "upset_lncrna",
                 text_panel("Fewer than two contrasts have significant lncRNAs."), 8, 6)
    }
  }, error = function(e) cat(sprintf("%s upset skipped: %s\n", sub, conditionMessage(e))))
}

# ---- updated wide table copy (consultation-friendly) ----
tryCatch({
  up <- ac
  sig_mat <- do.call(cbind, lapply(ct_ids, function(cid) sig_of(cid)))
  up$significant_contrasts <- apply(sig_mat, 1, function(r) paste(ct_ids[which(r)], collapse = ","))
  write_tsv(up, file.path(tables_dir, "lncrna_across_contrasts.tsv"))
}, error = function(e) cat("across table copy skipped:", conditionMessage(e), "\n"))

cat("contrast_summary_plots: done\n")
