#!/usr/bin/env Rscript
# Cross-contrast summary figures. Non-blocking: every panel is guarded so a
# missing package or an empty input never fails the task (the module also sets
# errorStrategy 'ignore'). Outputs both PDF (print) and PNG (slides).
#
#   contrast_summary_plots.R --overview contrasts_overview.tsv \
#       --across lncrna_across_contrasts.tsv --outdir .

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
overview_path <- getarg("--overview")
across_path   <- getarg("--across")
outdir        <- getarg("--outdir", ".")
fdr           <- as.numeric(getarg("--fdr", "0.05"))
dir.create(outdir, showWarnings = FALSE)

save_panel <- function(name, plot, width, height) {
  ggplot2::ggsave(file.path(outdir, paste0(name, ".png")), plot,
                  width = width, height = height, dpi = 150)
  ggplot2::ggsave(file.path(outdir, paste0(name, ".pdf")), plot,
                  width = width, height = height)
}

# ---- significant counts per contrast (up vs down) ----
tryCatch({
  ov <- read_tsv(overview_path, show_col_types = FALSE)
  if (nrow(ov) > 0 && all(c("n_sig_up", "n_sig_down") %in% colnames(ov))) {
    long <- ov %>%
      select(contrast, n_sig_up, n_sig_down) %>%
      pivot_longer(c(n_sig_up, n_sig_down), names_to = "dir", values_to = "n") %>%
      mutate(dir = recode(dir, n_sig_up = "up", n_sig_down = "down"),
             contrast = sub("__", " / ", contrast))
    p <- ggplot(long, aes(reorder(contrast, n, FUN = sum), n, fill = dir)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = c(up = "#b71c1c", down = "#1565c0")) +
      coord_flip() +
      labs(title = "Significant genes per contrast", x = NULL,
           y = "N significant (padj < FDR)", fill = NULL) +
      theme_minimal()
    save_panel("sig_counts_by_contrast", p, 9, max(4, 0.4 * nrow(ov) + 2))
  }
}, error = function(e) cat("sig_counts_by_contrast skipped:", conditionMessage(e), "\n"))

# ---- UpSet of lncRNAs significant across contrasts ----
tryCatch({
  ac <- read_tsv(across_path, show_col_types = FALSE)
  padj_cols <- grep("\\.padj$", colnames(ac), value = TRUE)
  if (nrow(ac) > 0 && length(padj_cols) >= 2) {
    memb <- as.data.frame(lapply(padj_cols, function(c) {
      as.integer(!is.na(ac[[c]]) & ac[[c]] < fdr)
    }))
    colnames(memb) <- sub("\\.padj$", "", padj_cols)
    # Drop contrasts with no significant lncRNA so an empty set does not crash upset()
    memb <- memb[, colSums(memb) > 0, drop = FALSE]
    if (ncol(memb) >= 2 && requireNamespace("UpSetR", quietly = TRUE)) {
      png(file.path(outdir, "upset_lncrna.png"), width = 1400, height = 800, res = 150)
      print(UpSetR::upset(memb, nsets = ncol(memb), order.by = "freq",
                          main.bar.color = "#37474f", sets.bar.color = "#c2185b"))
      dev.off()
      pdf(file.path(outdir, "upset_lncrna.pdf"), width = 12, height = 7)
      print(UpSetR::upset(memb, nsets = ncol(memb), order.by = "freq",
                          main.bar.color = "#37474f", sets.bar.color = "#c2185b"))
      dev.off()
    } else if (!requireNamespace("UpSetR", quietly = TRUE)) {
      cat("UpSetR not available; upset_lncrna skipped\n")
    }
  }
}, error = function(e) {
  # dev.off() may be pending if the plot failed mid-device
  try(grDevices::dev.off(), silent = TRUE)
  cat("upset_lncrna skipped:", conditionMessage(e), "\n")
})

cat("contrast_summary_plots: done\n")
