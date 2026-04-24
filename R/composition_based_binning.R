#' Composition based binning on metagenomic samples
#'
#' This function performs composition based binning on metagenomic samples,
#' directly from FASTA or FASTQ files, by short kmer analysis (k<8).
#' See \doi{10.1186/s12859-016-1186-3} for more details.
#'
#' @param ... Input sequences. Either character paths to FASTA/FASTQ files
#'   (uncompressed or gzip compressed), a \code{\link[Biostrings]{DNAStringSet}}
#'   / \code{\link[Biostrings]{QualityScaledDNAStringSet}}, or a
#'   \code{\link[ShortRead]{ShortReadQ}} object. Non-file inputs are staged to
#'   a temporary FASTA/FASTQ file for the Java backend.
#' @param kMerSizeCB kmer length for Composition based Binning.
#' @param numOfClustersCB Number of Clusters for Composition based Binning.
#' @param outputCB Output Composition based Binning Clusters
#'     files location and prefix.
#' @param keepQuality Keep fastq qualities on the output files.
#'     (will produce .fastq)
#' @param dryRun Don't write any output files.
#' @param gzip Gzip output files.
#' @param numOfThreads Number of threads to use. Defaults to
#'   \code{\link[BiocParallel]{bpworkers}()}.
#'
#' @return A \linkS4class{MetabinResult} object. Its \code{assignments} slot
#'   is a \code{\link[S4Vectors]{DataFrame}} with \code{numOfClustersCB + 2}
#'   columns:
#' \itemize{
#'     \item \code{read_id} : read identifier from fasta header
#'     \item \code{CB} : read was assigned to this CB cluster index
#'     \item \code{CB.n} : read to cluster CB.n distance
#' }
#' For backwards-compatible \code{data.frame} output use
#' \code{as.data.frame(result)}.
#' @export
#'
#' @examples
#' res <- composition_based_binning(
#'     system.file("extdata", "reads.metagenome.fasta.gz", package = "metabinR"),
#'     dryRun = TRUE, kMerSizeCB = 2
#' )
#' res
#' @author Anestis Gkanogiannis, \email{anestis@@gkanogiannis.com}
#' @references \url{https://github.com/gkanogiannis/metabinR}
composition_based_binning <- function(..., kMerSizeCB = 4,
                                      numOfClustersCB = 5,
                                      outputCB = "CB.cluster",
                                      keepQuality = FALSE, dryRun = FALSE,
                                      gzip = FALSE,
                                      numOfThreads = BiocParallel::bpworkers()) {
    inputs <- .resolve_inputs(list(...))

    params <- list(
        kMerSizeCB = kMerSizeCB, numOfClustersCB = numOfClustersCB,
        outputCB = outputCB,
        keepQuality = keepQuality, dryRun = dryRun, gzip = gzip,
        numOfThreads = numOfThreads
    )
    .check_binning_params(
        inputs = inputs, outputPrefix = outputCB,
        keepQuality = keepQuality, dryRun = dryRun, gzip = gzip,
        numOfThreads = numOfThreads, algo = "CB",
        extra = params
    )

    .call_bridge("CB", inputs = inputs, params = params)
}
