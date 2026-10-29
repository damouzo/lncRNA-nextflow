#!/usr/bin/env Rscript
# 2.1: DESeq2 must run on the tximport object so the average-transcript-length
# offset is kept. G1 has the same per-condition summed count but a very different
# average length, so the matrix path (counts only) sees almost no change while the
# offset path does.
#   apptainer exec $R_BIOC_IMG Rscript bin/tests/test_deseq2_offset.R

suppressPackageStartupMessages({ library(DESeq2); library(readr); library(tximport) })

args_all <- commandArgs(trailingOnly = FALSE)
self <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
self <- if (length(self) > 0) self[1] else "bin/tests/test_deseq2_offset.R"
script <- file.path(dirname(dirname(normalizePath(self))), "deseq2_contrasts.R")

fails <- 0
check <- function(cond, msg) {
  cat(sprintf("  %s %s\n", if (isTRUE(cond)) "ok  " else "FAIL", msg))
  if (!isTRUE(cond)) fails <<- fails + 1
}

work <- tempfile("off_"); dir.create(work)
old_wd <- getwd(); setwd(work)

samples <- c("A1", "A2", "A3", "B1", "B2", "B3")
cond <- rep(c("A", "B"), each = 3)
genes <- c("G1", "G2", "G3", "G4", "G5", "G6")

set.seed(1)
counts <- matrix(rnbinom(length(genes) * length(samples), mu = 200, size = 5),
                 nrow = length(genes), dimnames = list(genes, samples))
# G1: equal per-condition sums, values still vary across samples.
counts["G1", ] <- c(100, 110, 120, 100, 110, 120)
# G2: a stable housekeeping-like gene.
counts["G2", ] <- c(190, 200, 210, 190, 200, 210)

avg_len <- matrix(1000, nrow = length(genes), ncol = length(samples),
                  dimnames = list(genes, samples))
avg_len["G1", 1:3] <- 1800   # condition A uses the long isoform
avg_len["G1", 4:6] <- 520    # condition B uses the short isoform
abundance <- counts / avg_len
txi <- list(counts = counts, abundance = abundance, length = avg_len,
            countsFromAbundance = "no")

coldata <- data.frame(sample = samples, condition = cond, batch = "", norm_group = "G",
                      stringsAsFactors = FALSE)
write.csv(coldata, "coldata.csv", row.names = FALSE)
writeLines("A_vs_B\tA\tB\t", "contrasts.tsv")
saveRDS(counts, "counts.rds")

dds <- DESeqDataSetFromTximport(txi, coldata, ~condition)
saveRDS(dds, "dds.rds")

run <- function(extra) {
  system2("Rscript", c(script, "--group", "G", "--counts", "counts.rds",
                       extra, "--coldata", "coldata.csv",
                       "--contrasts", "contrasts.tsv",
                       "--design-formula", shQuote("~ condition"),
                       "--fdr", "0.05", "--outdir", "."))
}
off_flag <- function() {
  f <- "G_deseq2_results/G_length_offset.txt"
  if (!file.exists(f)) return(NA_character_)
  r <- read.delim(f, stringsAsFactors = FALSE)
  as.character(r$length_offset_used[1])
}

# Case 1: no --dds and no opt-in -> hard error.
st_no <- suppressWarnings(run(""))
check(st_no != 0, "no --dds and no --allow-no-offset fails (2.1)")

# Case 2: explicit opt-in -> runs and marks the offset as NOT used.
st_allow <- run("--allow-no-offset true")
de_mat <- as.data.frame(read_tsv("G_deseq2_results/A_vs_B_DE.tsv", show_col_types = FALSE))
check(st_allow == 0, "--allow-no-offset true runs")
check(identical(off_flag(), "FALSE"), "length_offset_used = FALSE without --dds")

# Case 3: with --dds -> runs and marks the offset as used.
st_off <- run("--dds dds.rds")
de_off <- as.data.frame(read_tsv("G_deseq2_results/A_vs_B_DE.tsv", show_col_types = FALSE))
check(st_off == 0, "--dds run exits 0")
check(identical(off_flag(), "TRUE"), "length_offset_used = TRUE with --dds")

G1m <- de_mat$log2FoldChange[match("G1", de_mat$gene_id)]
G1o <- de_off$log2FoldChange[match("G1", de_off$gene_id)]
cat(sprintf("G1 log2FC: matrix=%.4f offset=%.4f\n", G1m, G1o))
check(abs(G1m) < 0.5, "matrix path sees almost no change (summed counts equal)")
check(abs(G1o - G1m) > 0.5,
      "offset path accounts for the isoform-usage shift (LFC differs)")

setwd(old_wd)
cat("\n")
if (fails > 0) { cat(sprintf("FAILED: %d\n", fails)); quit(save = "no", status = 1) }
cat("ALL DESEQ2 OFFSET TESTS PASSED\n")
