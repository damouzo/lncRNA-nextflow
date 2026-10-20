// A2b: Strandedness input guard (fail-fast)
//
// Consumes the per-sample RSeQC outputs (previously computed but unused) and
// checks them against the declared library_strandedness contract BEFORE Salmon
// runs. Also reports the normalised fidelity E/(E+O) per sample.
//
// Fails when:
//   - samples disagree on paired/single-end library type
//   - the majority empirical orientation does not match library_strandedness
//     (declaring forward / unstranded on reverse data stops the run)
//   - strandedness_guard_strict=true and any sample is below the fidelity floor
// Warns (does not fail) when a sample is below strandedness_min_fidelity.
//
// Writes run_info/strandedness_summary.tsv. Never touches frozen_manifest.json.

process CHECK_STRANDEDNESS {

    input:
    path strandedness_files

    output:
    path "strandedness_summary.tsv", emit: summary

    script:
    def declared = Strandedness.mode(params.library_strandedness)
    def min_fid  = (params.strandedness_min_fidelity as double)
    def strict   = params.strandedness_guard_strict
    """
    #!/usr/bin/env Rscript

    DECLARED <- "${declared}"
    MIN_FID  <- as.numeric("${min_fid}")
    STRICT   <- ${strict ? 'TRUE' : 'FALSE'}

    PE_REV <- "1+-,1-+,2++,2--"
    PE_FWD <- "1++,1--,2+-,2-+"
    SE_REV <- "+-,-+"
    SE_FWD <- "++,--"

    files <- list.files(pattern = "\\\\.strandedness\\\\.txt\$", full.names = TRUE)
    if (length(files) == 0) stop("CHECK_STRANDEDNESS: no *.strandedness.txt files found")

    parse_one <- function(f) {
        lines <- readLines(f, warn = FALSE)
        paired <- any(grepl("PairEnd", lines, ignore.case = TRUE))
        single <- any(grepl("SingleEnd", lines, ignore.case = TRUE))
        data_type <- if (paired) "paired" else if (single) "single" else NA_character_

        frac_for <- function(pat) {
            hit <- grep(pat, lines, fixed = TRUE, value = TRUE)
            if (length(hit) == 0) return(NA_real_)
            as.numeric(trimws(sub(".*:", "", hit[1])))
        }
        failed_line <- grep("failed to determine", lines, fixed = TRUE, value = TRUE)
        failed <- if (length(failed_line) == 0) NA_real_ else
            as.numeric(trimws(sub(".*:", "", failed_line[1])))

        if (identical(data_type, "paired")) {
            rev <- frac_for(PE_REV); fwd <- frac_for(PE_FWD)
        } else {
            rev <- frac_for(SE_REV); fwd <- frac_for(SE_FWD)
        }
        list(sample = sub("\\\\.strandedness\\\\.txt\$", "", basename(f)),
             data_type = data_type, reverse = rev, forward = fwd, failed = failed)
    }

    parsed <- lapply(files, parse_one)
    bad <- vapply(parsed, function(p) is.na(p\$data_type) ||
                  is.na(p\$reverse) || is.na(p\$forward), logical(1))
    if (any(bad)) {
        stop(sprintf("CHECK_STRANDEDNESS: could not parse orientation for sample(s): %s",
                     paste(vapply(parsed[bad], function(p) p\$sample, character(1)), collapse = ", ")))
    }

    dtypes <- unique(vapply(parsed, function(p) p\$data_type, character(1)))
    if (length(dtypes) > 1) {
        stop(sprintf("CHECK_STRANDEDNESS: mixed library types across samples: %s",
                     paste(dtypes, collapse = ", ")))
    }

    rows <- lapply(parsed, function(p) {
        detected <- if (p\$reverse >= p\$forward) "reverse" else "forward"
        expected_orientation <- if (DECLARED == "unstranded") detected else DECLARED
        e <- if (expected_orientation == "reverse") p\$reverse else p\$forward
        o <- if (expected_orientation == "reverse") p\$forward else p\$reverse
        fidelity <- e / (e + o)
        data.frame(
            sample              = p\$sample,
            data_type           = p\$data_type,
            detected_orientation = detected,
            declared_strandedness = DECLARED,
            expected_fraction   = round(e, 4),
            contrary_fraction   = round(o, 4),
            failed_fraction     = round(p\$failed, 4),
            fidelity            = round(fidelity, 4),
            stringsAsFactors    = FALSE
        )
    })
    summ <- do.call(rbind, rows)
    summ <- summ[order(summ\$sample), ]
    write.table(summ, "strandedness_summary.tsv", sep = "\\t", quote = FALSE,
                row.names = FALSE)

    n_rev <- sum(summ\$detected_orientation == "reverse")
    n_fwd <- sum(summ\$detected_orientation == "forward")
    majority <- if (n_rev >= n_fwd) "reverse" else "forward"

    if (DECLARED == "unstranded") {
        stop(sprintf("CHECK_STRANDEDNESS: library_strandedness='unstranded' but data are %s-stranded (majority %s). Set library_strandedness accordingly.",
                     majority, majority))
    }
    if (majority != DECLARED) {
        stop(sprintf("CHECK_STRANDEDNESS: declared '%s' but empirical majority is '%s' (%d/%d samples). Fix library_strandedness or the samplesheet.",
                     DECLARED, majority, if (majority == 'reverse') n_rev else n_fwd, nrow(summ)))
    }

    low <- summ\$fidelity < MIN_FID
    if (any(low)) {
        msg <- sprintf("CHECK_STRANDEDNESS: %d sample(s) below fidelity floor %.2f: %s",
                       sum(low), MIN_FID, paste(summ\$sample[low], collapse = ", "))
        if (STRICT) stop(msg) else cat(sprintf("WARNING: %s\\n", msg))
    }
    cat(sprintf("Strandedness OK: %s (%d samples, median fidelity %.3f)\\n",
                DECLARED, nrow(summ), median(summ\$fidelity)))
    """
}
