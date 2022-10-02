#' Abundance based binning on metagenomic samples
#'
#' This function performs abundance based binning on metagenomic samples,
#' directly from FASTA or FASTQ files, by long kmer analysis (k>8).
#' See \doi{10.1186/s12859-016-1186-3} for more details.
#'
#' @param ... Input fasta/fastq files locations
#'     (uncompressed or gzip compressed).
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
#' @param numOfThreads Number of threads to use.
#'
#' @return A \code{\link[base]{data.frame}} of the binning assignments.
#'     Return value contains \code{numOfClustersAB + 2} columns.
#' \itemize{
#'     \item \code{read_id} : read identifier from fasta header
#'     \item \code{AB} : read was assigned to this AB cluster index
#'     \item \code{AB.n} : read to cluster AB.n distance
#' }
#' @export
#'
#' @examples
#' abundance_based_binning(
#'     system.file("extdata", "reads.metagenome.fasta.gz",package = "metabinR"),
#'     dryRun = TRUE, kMerSizeAB = 8
#' )
#' @author Anestis Gkanogiannis, \email{anestis@@gkanogiannis.com}
#' @references \url{https://github.com/gkanogiannis/metabinR}
#'

abundance_based_binning <- function(..., eMin = 1, eMax = 0, kMerSizeAB = 10,
                                    numOfClustersAB = 5, outputAB="AB.cluster",
                                    keepQuality = FALSE, dryRun = FALSE,
                                    gzip = FALSE, numOfThreads = 1) {
    ins <- unlist(list(...))

    abundance_based_binning_checkParams(ins = ins, eMin = eMin, eMax = eMax,
                                        kMerSizeAB = kMerSizeAB,
                                        numOfClustersAB = numOfClustersAB,
                                        outputAB = outputAB,
                                        keepQuality = keepQuality,
                                        dryRun = dryRun, gzip = gzip,
                                        numOfThreads = numOfThreads)

    metatarget <- rJava::.jnew(
        class="fr/cea/ig/metatarget/MTxAB",
        class.loader = .rJava.class.loader)
    cmd <- paste("--eMin", eMin, "--eMax", eMax, "--kMerSizeAB", kMerSizeAB,
                 "--numOfClustersAB", numOfClustersAB, "--outputAB", outputAB,
                 ifelse(keepQuality, "--quality", ""),
                 ifelse(dryRun, "--dry", ""),
                 ifelse(gzip, "--gzip", ""),
                 "--numOfThreads", numOfThreads,
                 "--input", paste(ins, collapse = " --input "), sep = " ")
    ret.str <- metatarget$go(rJava::.jarray(strsplit(cmd, "\\s+")[[1]]))

    return(utils::read.table(text = ret.str, header = TRUE))
}

abundance_based_binning_checkParams <- function(ins, eMin, eMax, kMerSizeAB,
                                                numOfClustersAB, outputAB,
                                                keepQuality, dryRun, gzip,
                                                numOfThreads) {
    if (length(ins)==0 || list(NULL) %in% ins) {
        stop("No input fasta/fastq files were provided.")
    }
    for (f in ins) {
        if(!methods::is(f, "character")){
            stop("Input file ", f, " does not exist.")
        }
        if(!file.exists(f)){
           stop("Input file ", f, " does not exist.")
        }
    }

    if (!is.numeric(eMin) ||
        (is.numeric(eMin) && eMin<1)) {
        stop("eMin parameter must be positive integer.")
    }

    if (!is.numeric(eMax) ||
        (is.numeric(eMax) && eMax<0)) {
        stop("eMax parameter must be integer >=0.")
    }

    if (!is.numeric(kMerSizeAB) || (is.numeric(kMerSizeAB) && kMerSizeAB<2)) {
        stop("kMerSizeAB parameter must be positive integer >1 .")
    }

    if (!is.numeric(numOfClustersAB) ||
        (is.numeric(numOfClustersAB) && numOfClustersAB<2)) {
        stop("numOfClustersAB parameter must be positive integer >1 .")
    }

    if ((!is.null(outputAB) && !methods::is(outputAB, "character")) ||
        (methods::is(outputAB, "character") && nchar(outputAB)==0)) {
        stop("outputAB must be a prefix location.")
    }

    if(!is.logical(keepQuality)) {
        stop("keepQuality parameter must be logical.")
    }

    if(!is.logical(dryRun)) {
        stop("dryRun parameter must be logical.")
    }

    if(!is.logical(gzip)) {
        stop("gzip parameter must be logical.")
    }

    if (!is.numeric(numOfThreads) ||
        (is.numeric(numOfThreads) && numOfThreads<1)) {
        stop("numOfThreads parameter must be positive integer.")
    }
}
