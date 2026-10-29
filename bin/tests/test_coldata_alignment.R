#!/usr/bin/env Rscript
# 0A.3: regression test for the H0 sample<->condition misalignment.
#
# Reproduces the NB4 case: the samplesheet lists the 5 conditions of a block in a
# fixed order (Scramble, DDX41sh1, DDX41sh2, DHX34sh1, DHX34sh2), but the quant
# directories are listed alphabetically by list.files, so colnames(counts) is a
# different order. Assigning rownames(coldata) <- colnames(counts) positionally
# then shifts every label; aligning by sample name does not.
#
# Runs bin/tximport_filter.R on synthetic quants and checks the published aligned
# coldata, plus the contrast that the shift would have broken.
#   apptainer exec $R_BIOC_IMG Rscript bin/tests/test_coldata_alignment.R

suppressPackageStartupMessages({ library(readr); library(DESeq2) })

args_all <- commandArgs(trailingOnly = FALSE)
self <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
self <- if (length(self) > 0) self[1] else "bin/tests/test_coldata_alignment.R"
script <- file.path(dirname(dirname(normalizePath(self))), "tximport_filter.R")
if (!file.exists(script)) stop("cannot locate bin/tximport_filter.R")

fails <- 0
check <- function(cond, msg) {
  cat(sprintf("  %s %s\n", if (isTRUE(cond)) "ok  " else "FAIL", msg))
  if (!isTRUE(cond)) fails <<- fails + 1
}

# Samplesheet order of one NB4 block (as written by the wet lab).
block <- c("Scramble", "DDX41sh1", "DDX41sh2", "DHX34sh1", "DHX34sh2")
samples <- paste0("NB4_Exp4_", block)
conditions <- paste0("NB4_", block)

work <- tempfile("cdal_"); dir.create(work)
old_wd <- getwd()
setwd(work)

# 300 genes; DDX41/DHX34 down in their knockdowns, everything else flat.
n_genes <- 300
genes <- sprintf("G%03d", 1:n_genes)
gids <- c("ENSG00000183258", "ENSG00000134815")
gene_ids <- c(gids, genes[-(1:2)])

tx <- data.frame(TXNAME = paste0(gene_ids, "_tx"), GENEID = gene_ids)
write_tsv(tx, "tx2gene.tsv")

expr <- function(gene_id, sample) {
  base <- 200
  if (gene_id == "ENSG00000183258") return(if (grepl("DDX41sh", sample)) 20 else 200)
  if (gene_id == "ENSG00000134815") return(if (grepl("DHX34sh", sample)) 20 else 200)
  base
}
for (s in samples) {
  qdir <- paste0(s, "_quant"); dir.create(qdir)
  num <- vapply(gene_ids, expr, numeric(1), sample = s)
  # per-sample and per-gene variation so VST has something to estimate
  num <- round(num * (1 + match(s, samples) * 0.07) * (1 + (seq_along(gene_ids) %% 7) * 0.02))
  write_tsv(data.frame(Name = tx$TXNAME, Length = 1000, EffectiveLength = 950,
                       TPM = num / 10, NumReads = num),
            file.path(qdir, "quant.sf"))
}

coldata <- data.frame(sample = samples, condition = conditions,
                      batch = "Exp4", norm_group = "NB4",
                      stringsAsFactors = FALSE)
write.csv(coldata, "coldata_NB4.csv", row.names = FALSE)

status <- system2("Rscript", c(script, "--group", "NB4", "--coldata", "coldata_NB4.csv",
                               "--tx2gene", "tx2gene.tsv", "--outdir", "."))
check(status == 0, "tximport_filter.R exits 0")
setwd(old_wd)

aligned_path <- file.path(work, "coldata_aligned_NB4.csv")
check(file.exists(aligned_path), "aligned coldata published")
if (file.exists(aligned_path)) {
  aligned <- read.csv(aligned_path, stringsAsFactors = FALSE)
  counts <- readRDS(file.path(work, "NB4_tximport_counts.rds"))

  check(identical(as.character(aligned$sample), colnames(counts)),
        "aligned coldata sample order == counts column order")
  check(identical(colnames(counts), sort(samples)),
        "counts columns are in alphabetical (list.files) order, not samplesheet order")

  # Each sample must keep the condition that belongs to its own name.
  check(identical(aligned$condition,
                  coldata$condition[match(aligned$sample, coldata$sample)]),
        "each sample keeps its own condition label")

  # The naive positional assignment must NOT be used: it reproduces the shift.
  naive <- setNames(coldata$condition, sort(samples))
  correct <- setNames(aligned$condition, aligned$sample)
  check(!identical(unname(naive[sort(samples)]), unname(correct[sort(samples)])),
        "positional assignment would shift the labels (bug reproduced)")
}

cat("\n")
if (fails > 0) { cat(sprintf("FAILED: %d\n", fails)); quit(save = "no", status = 1) }
cat("ALL COLDATA ALIGNMENT TESTS PASSED\n")
