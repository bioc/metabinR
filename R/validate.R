#' @importFrom checkmate assert_character assert_file_exists assert_int
#'   assert_count assert_flag assert_string
#' @importFrom cli cli_abort
NULL

# Coerce input paths from either `character(1+)` or a `list()` of scalars
# (the fallout of `abundance_based_binning(file1, file2, ...)` signatures).
.flatten_inputs <- function(ins) {
    if (length(ins) == 0L) {
        cli::cli_abort(
            "No input FASTA/FASTQ files were provided.",
            class = "metabinR_error_no_input"
        )
    }
    ins <- unlist(ins, use.names = FALSE)
    if (!is.character(ins)) {
        cli::cli_abort(
            "Inputs must be file paths given as {.cls character}.",
            class = "metabinR_error_input_type"
        )
    }
    ins
}

# Validate parameters common to all three binning algorithms, plus
# algorithm-specific parameters via the `extra` named list.
#
# `algo` is one of "AB", "CB", "ABxCB" and is used purely for error messages
# and to select which extra parameters are validated.
.check_binning_params <- function(inputs, outputPrefix, keepQuality, dryRun,
                                  gzip, numOfThreads, algo, extra = list()) {
    inputs <- .flatten_inputs(inputs)
    for (f in inputs) {
        checkmate::assert_file_exists(
            f, access = "r",
            .var.name = "input file"
        )
    }

    checkmate::assert_string(outputPrefix, min.chars = 1L,
                             .var.name = "output prefix")
    checkmate::assert_flag(keepQuality, .var.name = "keepQuality")
    checkmate::assert_flag(dryRun, .var.name = "dryRun")
    checkmate::assert_flag(gzip, .var.name = "gzip")
    checkmate::assert_count(numOfThreads, positive = TRUE,
                            .var.name = "numOfThreads")

    if (algo %in% c("AB", "ABxCB")) {
        checkmate::assert_count(extra$eMin, positive = TRUE,
                                .var.name = "eMin")
        checkmate::assert_int(extra$eMax, lower = 0L,
                              .var.name = "eMax")
        checkmate::assert_int(extra$kMerSizeAB, lower = 2L,
                              .var.name = "kMerSizeAB")
        checkmate::assert_int(extra$numOfClustersAB, lower = 2L,
                              .var.name = "numOfClustersAB")
    }
    if (algo %in% c("CB", "ABxCB")) {
        checkmate::assert_int(extra$kMerSizeCB, lower = 2L,
                              .var.name = "kMerSizeCB")
    }
    if (algo == "CB") {
        checkmate::assert_int(extra$numOfClustersCB, lower = 2L,
                              .var.name = "numOfClustersCB")
    }
    if (algo == "ABxCB") {
        checkmate::assert_count(extra$genomeSize, positive = TRUE,
                                .var.name = "genomeSize")
    }

    inputs
}
