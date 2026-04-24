#' Accessors for \linkS4class{MetabinResult}
#'
#' @param x A \linkS4class{MetabinResult} object.
#'
#' @return \code{assignments()} returns a \code{\link[S4Vectors]{DataFrame}}
#'   with the per-read cluster assignments and distances.
#'   \code{nClusters()} returns an integer scalar: the number of clusters
#'   inferred by the algorithm. \code{parameters()} returns the list of
#'   arguments passed to the binning function. \code{algorithm()} returns the
#'   algorithm tag (\code{"AB"}, \code{"CB"}, or \code{"ABxCB"}).
#'
#' @examples
#' res <- abundance_based_binning(
#'     system.file("extdata", "reads.metagenome.fasta.gz", package = "metabinR"),
#'     dryRun = TRUE, kMerSizeAB = 4, numOfClustersAB = 2
#' )
#' assignments(res)
#' nClusters(res)
#' algorithm(res)
#' parameters(res)$kMerSizeAB
#'
#' @name MetabinResult-accessors
#' @aliases assignments
#' @export
setGeneric("assignments", function(x) standardGeneric("assignments"))

#' @rdname MetabinResult-accessors
#' @export
setGeneric("nClusters", function(x) standardGeneric("nClusters"))

#' @rdname MetabinResult-accessors
#' @export
setGeneric("parameters", function(x) standardGeneric("parameters"))

#' @rdname MetabinResult-accessors
#' @export
setGeneric("algorithm", function(x) standardGeneric("algorithm"))
