#' Summarize the bins in a metabinR result
#'
#' Counts reads in each observed bin and summarizes their distance to the
#' assigned bin. Distances retain the scale of the binning algorithm; they are
#' not probabilities or comparable across algorithms.
#'
#' @param x A [MetabinResult-class] object.
#' @return A [S4Vectors::DataFrame] with `bin`, `n_reads`, `proportion`,
#'   `mean_distance`, and `median_distance`. Only observed bins are included.
#' @md
#' @export
#' @examples
#' res <- composition_based_binning(
#'     system.file("extdata", "reads.metagenome.fasta.gz", package = "metabinR"),
#'     dryRun = TRUE, kMerSizeCB = 2, numOfClustersCB = 2
#' )
#' bin_summary(res)
bin_summary <- function(x) {
    info <- .result_distance_columns(x)
    df <- info$assignments
    bins <- as.character(df[[info$algorithm]])
    levels <- unique(bins)
    n <- length(bins)

    counts <- tabulate(match(bins, levels), nbins = length(levels))
    means <- medians <- rep(NA_real_, length(levels))
    for (i in seq_along(levels)) {
        rows <- which(bins == levels[i])
        distance_col <- paste0(info$algorithm, ".", levels[i])
        if (!distance_col %in% info$distance_columns) {
            cli::cli_abort("No distance column for assigned bin {.val {levels[i]}}.")
        }
        distances <- df[[distance_col]][rows]
        distances <- distances[is.finite(distances)]
        if (length(distances)) {
            means[i] <- mean(distances)
            medians[i] <- stats::median(distances)
        }
    }

    S4Vectors::DataFrame(
        bin = levels,
        n_reads = counts,
        proportion = counts / n,
        mean_distance = means,
        median_distance = medians
    )
}

#' Find reads with close competing bins
#'
#' Returns reads for which the second-smallest finite distance is within
#' `margin` of the smallest distance. In hierarchical results, distances to
#' bins outside the read's parent abundance bin are `NA` and are ignored.
#' Reads with fewer than two finite distances are omitted. The margin is an
#' absolute difference on the algorithm's distance scale, not a probability.
#'
#' @param x A [MetabinResult-class] object.
#' @param margin Nonnegative maximum difference between the two smallest
#'   distances. Defaults to `0.05`.
#' @return A [S4Vectors::DataFrame] with `read_id`, `bin`, `best_distance`,
#'   `second_distance`, and `distance_margin`, in input order.
#' @md
#' @export
#' @examples
#' res <- composition_based_binning(
#'     system.file("extdata", "reads.metagenome.fasta.gz", package = "metabinR"),
#'     dryRun = TRUE, kMerSizeCB = 2, numOfClustersCB = 2
#' )
#' ambiguous_reads(res, margin = 0.05)
ambiguous_reads <- function(x, margin = 0.05) {
    if (!is.numeric(margin) || length(margin) != 1L ||
        !is.finite(margin) || margin < 0) {
        cli::cli_abort("{.arg margin} must be one finite nonnegative number.")
    }
    info <- .result_distance_columns(x)
    df <- info$assignments
    n <- nrow(df)
    best <- second <- rep(Inf, n)

    for (col in info$distance_columns) {
        distance <- df[[col]]
        valid <- is.finite(distance)
        better <- valid & distance < best
        second[better] <- best[better]
        best[better] <- distance[better]
        runner_up <- valid & !better & distance < second
        second[runner_up] <- distance[runner_up]
    }

    difference <- second - best
    selected <- which(is.finite(difference) & difference <= margin)
    S4Vectors::DataFrame(
        read_id = df[["read_id"]][selected],
        bin = as.character(df[[info$algorithm]][selected]),
        best_distance = best[selected],
        second_distance = second[selected],
        distance_margin = difference[selected]
    )
}

.result_distance_columns <- function(x) {
    if (!methods::is(x, "MetabinResult")) {
        cli::cli_abort("{.arg x} must be a {.cls MetabinResult} object.")
    }
    df <- assignments(x)
    algo <- algorithm(x)
    if (is.na(algo) || !all(c("read_id", algo) %in% names(df))) {
        cli::cli_abort("Result has no usable read and bin assignments.")
    }
    if (anyNA(df[[algo]])) {
        cli::cli_abort("Result has missing bin assignments.")
    }
    cols <- grep(paste0("^", algo, "\\.[0-9]+$"), names(df), value = TRUE)
    if (!length(cols) ||
        any(!vapply(cols, function(col) is.numeric(df[[col]]), logical(1)))) {
        cli::cli_abort("Result has no usable numeric bin distances.")
    }
    list(assignments = df, algorithm = algo, distance_columns = cols)
}
