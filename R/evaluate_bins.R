#' Evaluate read-level bin assignments against known origins
#'
#' Joins assignments to a truth table by read identifier, then calculates
#' per-bin purity, per-origin best-bin recovery, and the adjusted Rand index
#' (ARI). These are read-count metrics. They do not measure genome or
#' metagenome-assembled genome completeness, and they are not CAMI/AMBER
#' base-pair-weighted scores. Bin labels are treated as arbitrary identifiers.
#'
#' Every result read must have exactly one matching truth row. Extra truth
#' rows are ignored. Duplicate or missing identifiers are errors. If several
#' truth groups or bins tie for the maximum count, the first encountered is
#' reported as dominant.
#'
#' @param x A [MetabinResult-class] object.
#' @param truth A `data.frame` or [S4Vectors::DataFrame] containing read IDs
#'   and known origins or classes.
#' @param id Name of the read-ID column in `truth`. The result's `read_id`
#'   column is matched to this column.
#' @param label Name of the truth-label column in `truth`.
#' @return A list with `confusion` (bins as rows, truth groups as columns),
#'   `per_bin` (read count, dominant truth group, purity), `per_origin`
#'   (read count, dominant bin, best-bin recovery), and `overall` (read count,
#'   numbers of bins and truth groups, weighted purity, weighted recovery,
#'   and ARI). Purity is the dominant truth count divided by bin size;
#'   best-bin recovery is the largest count for an origin in any one bin
#'   divided by that origin's read count. Both weighted measures sum the
#'   relevant dominant counts and divide by the total number of evaluated
#'   reads. ARI is `NA` for fewer than two reads.
#' @md
#' @export
#' @examples
#' x <- composition_based_binning(
#'     system.file("extdata", "reads.metagenome.fasta.gz", package = "metabinR"),
#'     dryRun = TRUE, kMerSizeCB = 2, numOfClustersCB = 2
#' )
#' truth <- read.delim(system.file(
#'     "extdata", "reads_mapping.tsv.gz", package = "metabinR"
#' ))
#' evaluation <- evaluate_bins(x, truth, id = "anonymous_read_id",
#'                             label = "genome_id")
#' evaluation$overall
evaluate_bins <- function(x, truth, id = "read_id", label = "genome_id") {
    if (!methods::is(x, "MetabinResult")) {
        cli::cli_abort("{.arg x} must be a {.cls MetabinResult} object.")
    }
    if (!is.data.frame(truth) && !methods::is(truth, "DataFrame")) {
        cli::cli_abort("{.arg truth} must be a data.frame or DataFrame.")
    }
    if (!is.character(id) || length(id) != 1L || is.na(id) || !nzchar(id) ||
        !is.character(label) || length(label) != 1L || is.na(label) ||
        !nzchar(label)) {
        cli::cli_abort("{.arg id} and {.arg label} must be column names.")
    }
    if (!all(c(id, label) %in% names(truth))) {
        cli::cli_abort("Truth table lacks the requested ID or label column.")
    }

    df <- assignments(x)
    algo <- algorithm(x)
    if (is.na(algo) || !all(c("read_id", algo) %in% names(df))) {
        cli::cli_abort("Result has no usable read and bin assignments.")
    }
    read_ids <- as.character(df[["read_id"]])
    truth_ids <- as.character(truth[[id]])
    if (!length(read_ids) || anyNA(read_ids) || any(!nzchar(read_ids)) ||
        anyDuplicated(read_ids)) {
        cli::cli_abort("Result read IDs must be nonempty, unique, and nonmissing.")
    }
    if (anyNA(truth_ids) || any(!nzchar(truth_ids)) ||
        anyDuplicated(truth_ids)) {
        cli::cli_abort("Truth read IDs must be nonempty, unique, and nonmissing.")
    }
    matched <- match(read_ids, truth_ids)
    if (anyNA(matched)) {
        cli::cli_abort("Truth table is missing {sum(is.na(matched))} result read(s).")
    }
    bins <- as.character(df[[algo]])
    origins <- as.character(truth[[label]][matched])
    if (anyNA(bins) || any(!nzchar(bins)) ||
        anyNA(origins) || any(!nzchar(origins))) {
        cli::cli_abort("Bin and truth labels must be nonempty and nonmissing.")
    }

    bin_levels <- unique(bins)
    origin_levels <- unique(origins)
    confusion <- table(
        bin = factor(bins, levels = bin_levels),
        origin = factor(origins, levels = origin_levels)
    )
    bin_counts <- rowSums(confusion)
    origin_counts <- colSums(confusion)
    bin_dominant_index <- max.col(confusion, ties.method = "first")
    origin_dominant_index <- max.col(t(confusion), ties.method = "first")
    bin_dominant_counts <- confusion[
        cbind(seq_along(bin_levels), bin_dominant_index)
    ]
    origin_dominant_counts <- confusion[
        cbind(origin_dominant_index, seq_along(origin_levels))
    ]

    per_bin <- S4Vectors::DataFrame(
        bin = bin_levels,
        n_reads = as.integer(bin_counts),
        dominant_origin = origin_levels[bin_dominant_index],
        dominant_reads = as.integer(bin_dominant_counts),
        purity = as.numeric(bin_dominant_counts / bin_counts)
    )
    per_origin <- S4Vectors::DataFrame(
        origin = origin_levels,
        n_reads = as.integer(origin_counts),
        dominant_bin = bin_levels[origin_dominant_index],
        recovered_reads = as.integer(origin_dominant_counts),
        best_bin_recovery = as.numeric(origin_dominant_counts / origin_counts)
    )
    n <- length(read_ids)
    overall <- S4Vectors::DataFrame(
        n_reads = n,
        n_bins = length(bin_levels),
        n_origins = length(origin_levels),
        weighted_purity = sum(bin_dominant_counts) / n,
        weighted_recovery = sum(origin_dominant_counts) / n,
        adjusted_rand_index = .adjusted_rand_index(confusion)
    )
    list(confusion = confusion, per_bin = per_bin,
         per_origin = per_origin, overall = overall)
}

.adjusted_rand_index <- function(confusion) {
    n <- sum(confusion)
    if (n < 2L) return(NA_real_)
    pairs <- function(counts) {
        counts <- as.numeric(counts)
        sum(counts * (counts - 1) / 2)
    }
    total_pairs <- n * (n - 1) / 2
    within_cells <- pairs(confusion)
    within_bins <- pairs(rowSums(confusion))
    within_origins <- pairs(colSums(confusion))
    expected <- within_bins * within_origins / total_pairs
    maximum <- (within_bins + within_origins) / 2
    denominator <- maximum - expected
    if (denominator == 0) return(1)
    (within_cells - expected) / denominator
}
