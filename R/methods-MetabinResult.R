#' @rdname MetabinResult-accessors
#' @export
setMethod("assignments", "MetabinResult", function(x) x@assignments)

#' @rdname MetabinResult-accessors
#' @export
setMethod("parameters", "MetabinResult", function(x) x@parameters)

#' @rdname MetabinResult-accessors
#' @export
setMethod("algorithm", "MetabinResult", function(x) x@algorithm)

#' @rdname MetabinResult-accessors
#' @export
setMethod("nClusters", "MetabinResult", function(x) {
    df <- x@assignments
    col <- x@algorithm
    if (length(col) != 1L || is.na(col) || !col %in% colnames(df)) {
        return(0L)
    }
    length(unique(df[[col]]))
})

#' Convert a \linkS4class{MetabinResult} to a \code{data.frame}
#'
#' Preserves backwards compatibility with the \code{data.frame} return of
#' metabinR <= 1.x.
#'
#' @param x A \linkS4class{MetabinResult}.
#' @param row.names,optional,... Unused; retained for S3 signature compatibility.
#' @return A base \code{data.frame} of the assignments.
#' @export
setMethod(
    "as.data.frame", "MetabinResult",
    function(x, row.names = NULL, optional = FALSE, ...) {
        as.data.frame(x@assignments, row.names = row.names, optional = optional)
    }
)

#' @rdname MetabinResult-class
#' @param object A \linkS4class{MetabinResult}.
#' @export
setMethod("show", "MetabinResult", function(object) {
    n <- nrow(object@assignments)
    k <- nClusters(object)
    cat("MetabinResult (", object@algorithm, ")\n", sep = "")
    cat("  reads:     ", format(n, big.mark = ","), "\n", sep = "")
    cat("  clusters:  ", k, "\n", sep = "")
    cat("  inputs:    ", length(object@inputs), " file(s)\n", sep = "")
    if (length(object@inputs) > 0L && length(object@inputs) <= 3L) {
        for (f in object@inputs) cat("    - ", f, "\n", sep = "")
    }
    invisible(NULL)
})
