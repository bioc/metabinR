#' Hierarchical (ABxCB) binning on metagenomic samples
#'
#' This function performs hierarchical binning on metagenomic samples,
#' directly from FASTA or FASTQ files.
#' First it analyzes sequences by long kmer analysis (k>8),
#' as in \code{\link[metabinR]{abundance_based_binning}}.
#' Then for each AB bin, it guesses the number of composition bins in it and
#' performs composition based binning by short kmer analysis (k<8),
#' as in \code{\link[metabinR]{composition_based_binning}}.
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
#' @param kMerSizeCB kmer length for Composition based Binning.
#' @param genomeSize Average genome size of taxa in the metagenome data.
#' @param numOfClustersAB Number of Clusters for Abundance based Binning.
#' @param outputC Output Hierarchical Binning (ABxCB) Clusters
#'     files location and prefix.
#' @param keepQuality Keep fastq qualities on the output files.
#'     (will produce .fastq)
#' @param dryRun Don't write any output files.
#' @param gzip Gzip output files.
#' @param numOfThreads Number of threads to use. Defaults to
#'   \code{\link[BiocParallel]{bpworkers}()}.
#'
#' @return A \linkS4class{MetabinResult} object. Its \code{assignments} slot
#'   is a \code{\link[S4Vectors]{DataFrame}}:
#' \itemize{
#'     \item \code{read_id} : read identifier from fasta header
#'     \item \code{ABxCB} : read was assigned to this ABxCB cluster index
#'     \item \code{ABxCB.n} : read to cluster ABxCB.n distance
#' }
#' For backwards-compatible \code{data.frame} output use
#' \code{as.data.frame(result)}.
#' @export
#'
#' @examples
#' res <- hierarchical_binning(
#'     system.file("extdata", "reads.metagenome.fasta.gz", package = "metabinR"),
#'     dryRun = TRUE, kMerSizeAB = 4, kMerSizeCB = 2
#' )
#' res
#' @author Anestis Gkanogiannis, \email{anestis@@gkanogiannis.com}
#' @references \url{https://github.com/gkanogiannis/metabinR}
hierarchical_binning <- function(..., eMin = 1, eMax = 0, kMerSizeAB = 10,
                                 kMerSizeCB = 4, genomeSize = 3000000,
                                 numOfClustersAB = 3,
                                 outputC = "ABxCB.cluster",
                                 keepQuality = FALSE, dryRun = FALSE,
                                 gzip = FALSE,
                                 numOfThreads = BiocParallel::bpworkers()) {
    inputs <- .resolve_inputs(list(...))

    params <- list(
        eMin = eMin, eMax = eMax,
        kMerSizeAB = kMerSizeAB, kMerSizeCB = kMerSizeCB,
        genomeSize = genomeSize,
        numOfClustersAB = numOfClustersAB,
        readLength = 0L,
        outputC = outputC,
        keepQuality = keepQuality, dryRun = dryRun, gzip = gzip,
        numOfThreads = numOfThreads
    )
    .check_binning_params(
        inputs = inputs, outputPrefix = outputC,
        keepQuality = keepQuality, dryRun = dryRun, gzip = gzip,
        numOfThreads = numOfThreads, algo = "ABxCB",
        extra = params
    )

    .call_bridge("ABxCB", inputs = inputs, params = params)
}
