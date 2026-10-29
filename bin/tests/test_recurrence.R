#!/usr/bin/env Rscript
# 3.1: recurrence CPM denominator must be the total fragments of the sample
# (featureCounts .summary, all categories), never the candidate-only column sums.
#
# Two overlapping isoforms plus a large Unassigned_NoFeatures block make the two
# denominators disagree: txA has CPM 2 with the real denominator (fails a 5 CPM
# cutoff) but a huge CPM with the candidate sum. The test asserts the real
# behaviour and, by dropping the summaries, shows the naive one differing.
#   apptainer exec $R_BIOC_IMG Rscript bin/tests/test_recurrence.R

suppressPackageStartupMessages({ library(readr) })

args_all <- commandArgs(trailingOnly = FALSE)
self <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
self <- if (length(self) > 0) self[1] else "bin/tests/test_recurrence.R"
script <- file.path(dirname(dirname(normalizePath(self))), "expression_recurrence.R")

fails <- 0
check <- function(cond, msg) {
  cat(sprintf("  %s %s\n", if (isTRUE(cond)) "ok  " else "FAIL", msg))
  if (!isTRUE(cond)) fails <<- fails + 1
}

work <- tempfile("rec_"); dir.create(work)
old_wd <- getwd(); setwd(work)

writeLines(c(
  'chr1\tStringTie\ttranscript\t100\t200\t.\t+\t.\tgene_id "gA"; transcript_id "txA"; class_code "u";',
  'chr1\tStringTie\ttranscript\t150\t400\t.\t+\t.\tgene_id "gB"; transcript_id "txB"; class_code "u";'
), "merged.gtf")
writeLines("transcript_id\tgene_id\ntxA\tgA\ntxB\tgB", "catalog.tsv")

# Candidate counts: txA = 2, txB = 100 in each sample. A large unassigned block
# keeps the total fragments near 1e6, so txA's real CPM is 2 (fails a 5 cutoff).
mk_fc <- function(sample) {
  writeLines(c(
    "# Program:featureCounts v2.0.1",
    sprintf("Geneid\tChr\tStart\tEnd\tStrand\tLength\t%s.bam", sample),
    "txA\tchr1\t100\t200\t+\t101\t2",
    "txB\tchr1\t150\t400\t+\t251\t100"
  ), sprintf("%s.featureCounts.tsv", sample))
  # Total fragments of the sample = sum of every category.
  writeLines(c(
    sprintf("Status\t%s.bam", sample),
    "Assigned\t102",
    "Unassigned_Unmapped\t0",
    sprintf("Unassigned_NoFeatures\t%d", 1000000 - 102)
  ), sprintf("%s.featureCounts.tsv.summary", sample))
}
for (s in c("S1", "S2")) {
  mk_fc(s)
}

run <- function(...) system2("Rscript", c(script, ..., "--outdir", "."))
st <- suppressWarnings(run("--catalog", "catalog.tsv", "--merged-gtf", "merged.gtf",
                           "--fc-dir", ".", "--min-reps", "2", "--cpm-threshold", "5"))
check(st == 0 && file.exists("recurrence_filtered.tsv"), "expression_recurrence.R exits 0")

rec <- read_tsv("recurrence_filtered.tsv", show_col_types = FALSE)
pr <- setNames(rec$passes_recurrence, rec$transcript_id)
check(isFALSE(pr[["txA"]]), "txA fails recurrence with the total-fragment denominator (CPM=2)")
check(isTRUE(pr[["txB"]]), "txB passes recurrence (CPM=100)")

stats <- readLines("recurrence_stats.txt")
check(any(grepl("cpm_denominator\tper_sample_total_fragments", stats)),
      "effective denominator registered in recurrence_stats.txt")
check(any(grepl("cpm_threshold\t5.0", stats)), "effective CPM threshold registered")

# Without the summaries the naive candidate-sum denominator flips txA to pass.
file.remove(list.files(pattern = "\\.summary$", full.names = TRUE))
st2 <- suppressWarnings(run("--catalog", "catalog.tsv", "--merged-gtf", "merged.gtf",
                            "--fc-dir", ".", "--min-reps", "2", "--cpm-threshold", "5"))
rec2 <- read_tsv("recurrence_filtered.tsv", show_col_types = FALSE)
pr2 <- setNames(rec2$passes_recurrence, rec2$transcript_id)
check(st2 == 0 && isTRUE(pr2[["txA"]]),
      "fallback (no summary) inflates txA CPM -> the denominator really matters")

setwd(old_wd)
cat("\n")
if (fails > 0) { cat(sprintf("FAILED: %d\n", fails)); quit(save = "no", status = 1) }
cat("ALL RECURRENCE TESTS PASSED\n")
