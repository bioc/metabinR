#' MetabinResult: binning result container
#'
#' An S4 class returned by [abundance_based_binning()],
#' [composition_based_binning()] and [hierarchical_binning()].
#'
#' @slot assignments A \code{\link[S4Vectors]{DataFrame}} of cluster
#'     assignments. The first column is \code{read_id}; subsequent columns
#'     are algorithm-specific (see the corresponding binning function).
#' @slot parameters Named list of the parameters passed to the algorithm.
#' @slot inputs Character vector of input file paths that were processed.
#' @slot algorithm Character scalar: one of \code{"AB"}, \code{"CB"},
#'     \code{"ABxCB"}.
#'
#' @return Objects of this class are returned by the binning functions.
#'   Use \code{\link{assignments}}, \code{\link{nClusters}},
#'   \code{\link{parameters}}, \code{\link{algorithm}}, or
#'   \code{as.data.frame()} to access the results.
#' @examples
#' res <- abundance_based_binning(
#'     system.file("extdata", "reads.metagenome.fasta.gz", package = "metabinR"),
#'     dryRun = TRUE, kMerSizeAB = 4, numOfClustersAB = 2
#' )
#' res
#' is(res, "MetabinResult")
#'
#' @name MetabinResult-class
#' @aliases MetabinResult
#' @exportClass MetabinResult
setClass(
    "MetabinResult",
    representation(
        assignments = "DataFrame",
        parameters = "list",
        inputs = "character",
        algorithm = "character"
    ),
    prototype(
        assignments = S4Vectors::DataFrame(),
        parameters = list(),
        inputs = character(0),
        algorithm = NA_character_
    )
)

setValidity("MetabinResult", function(object) {
    errs <- character(0)
    if (length(object@algorithm) != 1L ||
        !object@algorithm %in% c("AB", "CB", "ABxCB", NA_character_)) {
        errs <- c(errs, "`algorithm` must be one of 'AB', 'CB', 'ABxCB'.")
    }
    if (length(errs) == 0L) TRUE else errs
})
