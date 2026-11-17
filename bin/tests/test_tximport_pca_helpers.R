#!/usr/bin/env Rscript
# Helpers for the per-group PCA panels in bin/tximport_filter.R (palette + plot
# builder). Parses the real script and evaluates only the helper definitions, so
# it does not need DESeq2/limma and runs cheaply in the r_bioc image:
#   apptainer exec $IMG Rscript bin/tests/test_tximport_pca_helpers.R

suppressPackageStartupMessages({library(ggplot2); library(readr)})

args_all <- commandArgs(trailingOnly = FALSE)
self <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
self <- if (length(self) > 0) self[1] else "bin/tests/test_tximport_pca_helpers.R"
script <- file.path(dirname(dirname(normalizePath(self))), "tximport_filter.R")
if (!file.exists(script)) stop("cannot locate bin/tximport_filter.R next to the test")

fails <- 0
check <- function(cond, msg) {
  cat(sprintf("  %s %s\n", if (isTRUE(cond)) "ok  " else "FAIL", msg))
  if (!isTRUE(cond)) fails <<- fails + 1
}

# Evaluate only the helper assignments (palette + plot builder), not the pipeline.
want <- c("CORONA", "BREWER_MAX", "big_palette", "pca_scatter", "write_pca", "has_levels")
env <- new.env()
for (e in parse(script)) {
  if (is.call(e) && length(e) == 3 && identical(e[[1]], as.name("<-")) &&
      is.symbol(e[[2]]) && as.character(e[[2]]) %in% want) {
    eval(e, env)
  }
}
check(all(want %in% ls(env)),
      sprintf("helpers parsed from the script (%s)", paste(sort(ls(env)), collapse = ", ")))

# Palette: >= n unique hex colours for a full cohort, from corona+Set1+Dark2.
for (n in c(3, 10, 38, 60)) {
  p <- env$big_palette(n, "corona+Set1+Dark2")
  check(length(p) == n && length(unique(p)) == n,
        sprintf("big_palette(%d) -> %d unique colours", n, length(unique(p))))
  check(all(grepl("^#[0-9A-Fa-f]{6}$", p)), sprintf("big_palette(%d) all hex", n))
}
check(identical(env$CORONA[1], "#1f77b4") && length(env$CORONA) == 30,
      "corona palette (BuenColors) has its 30 colours")

# has_levels drops single-level / empty columns (batch is '' in primary groups).
check(!env$has_levels(c("", "", "")), "has_levels: empty column is single level")
check(!env$has_levels(c("main", "main")), "has_levels: constant column is single level")
check(env$has_levels(c("main", "L4")), "has_levels: two values -> TRUE")

# pca_scatter builds a ggplot with colour and shape kept as separate aes.
df <- data.frame(PC1 = rnorm(6), PC2 = rnorm(6), name = letters[1:6], check.names = FALSE)
df$.color_group <- factor(c("a", "a", "b", "b", "b", "a"))
df$.shape_group <- factor(c("x", "y", "x", "y", "x", "y"))
g <- env$pca_scatter(df, c(40L, 20L),
                     ggplot2::aes(x = PC1, y = PC2, color = .color_group, shape = .shape_group),
                     "t", "s", "Condition", "Batch")
check(inherits(g, "ggplot"), "pca_scatter returns a ggplot")
built <- ggplot2::ggplot_build(g)
check(length(built$data) >= 1, "pca_scatter renders (ggplot_build)")
check(identical(built$plot$labels$colour, "Condition") &&
        identical(built$plot$labels$shape, "Batch"),
      "pca_scatter keeps shape and colour labels separate")

g2 <- env$pca_scatter(df, c(40L, 20L),
                      ggplot2::aes(x = PC1, y = PC2, color = .color_group),
                      "t", "s", "Condition")
check(inherits(g2, "ggplot"), "pca_scatter colour-only returns a ggplot")
invisible(ggplot2::ggplot_build(g2))

cat("\n")
if (fails > 0) { cat(sprintf("FAILED: %d\n", fails)); quit(save = "no", status = 1) }
cat("ALL PCA HELPER TESTS PASSED\n")
