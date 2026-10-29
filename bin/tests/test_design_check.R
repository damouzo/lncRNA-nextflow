#!/usr/bin/env Rscript
# F3.2 test: design_check.R on synthetic metadata.
#   apptainer exec $R_BIOC_IMG Rscript bin/tests/test_design_check.R

suppressPackageStartupMessages(library(readr))

args_all <- commandArgs(trailingOnly = FALSE)
self <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
self <- if (length(self) > 0) self[1] else "bin/tests/test_design_check.R"
script <- file.path(dirname(dirname(normalizePath(self))), "design_check.R")
if (!file.exists(script)) stop("cannot locate bin/design_check.R")

fails <- 0
check <- function(cond, msg) {
  cat(sprintf("  %s %s\n", if (isTRUE(cond)) "ok  " else "FAIL", msg))
  if (!isTRUE(cond)) fails <<- fails + 1
}

run <- function(coldata, contrasts, mode) {
  d <- tempfile("dc_"); dir.create(d)
  write.csv(coldata, file.path(d, "coldata.csv"), row.names = FALSE)
  writeLines(contrasts, file.path(d, "contrasts.tsv"))
  st <- system2("Rscript", c(script, "--group", "G", "--coldata", file.path(d, "coldata.csv"),
                             "--contrasts", file.path(d, "contrasts.tsv"),
                             "--design-formula", shQuote("~ batch + condition"),
                             "--mode", mode, "--outdir", d))
  list(status = st, tsv = file.path(d, "design_check.tsv"))
}

# Clean design: batch balances condition -> full rank
clean <- data.frame(
  sample = paste0("S", 1:6),
  condition = c("A", "A", "A", "B", "B", "B"),
  batch = c("b1", "b1", "b2", "b1", "b2", "b2"),
  stringsAsFactors = FALSE)
r <- run(clean, "A_vs_B\tA\tB\t", "fail")
check(r$status == 0, "clean design passes")
if (file.exists(r$tsv)) {
  dc <- read_tsv(r$tsv, show_col_types = FALSE)
  check(all(dc$status == "ok"), "clean design status ok")
}

# Confounded: every A in b1, every B in b2 -> ~ batch + condition rank-deficient
conf <- data.frame(
  sample = paste0("S", 1:6),
  condition = c("A", "A", "A", "B", "B", "B"),
  batch = c("b1", "b1", "b1", "b2", "b2", "b2"),
  stringsAsFactors = FALSE)
r <- run(conf, "A_vs_B\tA\tB\t", "fail")
check(r$status != 0, "confounded design fails with mode=fail")

r <- run(conf, "A_vs_B\tA\tB\t", "warn")
check(r$status == 0, "confounded design continues with mode=warn")
if (file.exists(r$tsv)) {
  dc <- read_tsv(r$tsv, show_col_types = FALSE)
  check(any(dc$status == "rank_deficient"), "confounded design recorded as rank_deficient")
}

cat("\n")
if (fails > 0) { cat(sprintf("FAILED: %d\n", fails)); quit(save = "no", status = 1) }
cat("ALL DESIGN CHECK TESTS PASSED\n")
