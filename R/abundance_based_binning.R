#' Abundance based binning on metagenomic samples
#'
#' This function performs abundance based binning on metagenomic samples,
#' directly from FASTA or FASTQ files, by long kmer analysis (k>8).
#' See \doi{10.1186/s12859-016-1186-3} for more details.
#'
#' @param ... Input sequences. Either character paths to FASTA/FASTQ files
#'   (uncompressed or gzip compressed), a \code{\link[Biostrings]{DNAStringSet}}
#'   / \code{\link[Biostrings]{QualityScaledDNAStringSet}}, or a
#'   \code{\link[ShortRead]{ShortReadQ}} object. Non-file inputs are staged to
#'   a temporary FASTA/FASTQ file for the Java backend.
#' @param eMin Exclude kmers of less or equal count.
#' @param eMax Exclude kmers of more or equal count.
#' @param kMerSizeAB kmer length for Abundance based Binning.
#' @param numOfClustersAB Number of Clusters for Abundance based Binning.
#' @param outputAB Output Abundance based Binning Clusters
#'     files location and prefix.
#' @param keepQuality Keep fastq qualities on the output files.
#'     (will produce .fastq)
#' @param dryRun Don't write any output files.
#' @param gzip Gzip output files.
#' @param numOfThreads Number of threads to use. Defaults to
#'   \code{\link[BiocParallel]{bpworkers}()}.
#'
#' @return A \linkS4class{MetabinResult} object. Its \code{assignments} slot
#'   is a \code{\link[S4Vectors]{DataFrame}} with \code{numOfClustersAB + 2}
#'   columns:
#' \itemize{
#'     \item \code{read_id} : read identifier from fasta header
#'     \item \code{AB} : read was assigned to this AB cluster index
#'     \item \code{AB.n} : read to cluster AB.n distance
#' }
#' For backwards-compatible \code{data.frame} output use
#' \code{as.data.frame(result)}.
#' @export
#'
#' @examples
#' res <- abundance_based_binning(
#'     system.file("extdata", "reads.metagenome.fasta.gz", package = "metabinR"),
#'     dryRun = TRUE, kMerSizeAB = 8
#' )
#' res
#' head(as.data.frame(res))
#' @author Anestis Gkanogiannis, \email{anestis@@gkanogiannis.com}
#' @references \url{https://github.com/gkanogiannis/metabinR}
abundance_based_binning <- function(..., eMin = 1, eMax = 0, kMerSizeAB = 10,
                                    numOfClustersAB = 3,
                                    outputAB = "AB.cluster",
                                    keepQuality = FALSE, dryRun = FALSE,
                                    gzip = FALSE,
                                    numOfThreads = BiocParallel::bpworkers()) {
    inputs <- .resolve_inputs(list(...))

    params <- list(
        eMin = eMin, eMax = eMax,
        kMerSizeAB = kMerSizeAB, numOfClustersAB = numOfClustersAB,
        outputAB = outputAB,
        keepQuality = keepQuality, dryRun = dryRun, gzip = gzip,
        numOfThreads = numOfThreads
    )
    .check_binning_params(
        inputs = inputs, outputPrefix = outputAB,
        keepQuality = keepQuality, dryRun = dryRun, gzip = gzip,
        numOfThreads = numOfThreads, algo = "AB",
        extra = params
    )

    .call_bridge("AB", inputs = inputs, params = params)
}
