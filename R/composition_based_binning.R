#' Composition based binning on metagenomic samples
#'
#' This function performs composition based binning on metagenomic samples,
#' directly from FASTA or FASTQ files, by short kmer analysis (k<8).
#' See \doi{10.1186/s12859-016-1186-3} for more details.
#'
#' @param ... Input fasta/fastq files locations
#'     (uncompressed or gzip compressed).
#' @param kMerSizeCB kmer length for Composition based Binning.
#' @param numOfClustersCB Number of Clusters for Composition based Binning.
#' @param outputCB Output Composition based Binning Clusters
#'     files location and prefix.
#' @param keepQuality Keep fastq qualities on the output files.
#'     (will produce .fastq)
#' @param dryRun Don't write any output files.
#' @param gzip Gzip output files.
#' @param numOfThreads Number of threads to use.
#'
#' @return A \code{\link[base]{data.frame}} of the binning assignments.
#'     Return value contains \code{numOfClustersCB + 2} columns.
#' \itemize{
#'     \item \code{read_id} : read identifier from fasta header
#'     \item \code{CB} : read was assigned to this CB cluster index
#'     \item \code{CB.n} : read to cluster CB.n distance
#' }
#' @export
#'
#' @examples
#' composition_based_binning(
#'     system.file("extdata", "reads.metagenome.fasta.gz",package = "metabinR"),
#'     dryRun = TRUE, kMerSizeCB = 2
#' )
#' @author Anestis Gkanogiannis, \email{anestis@@gkanogiannis.com}
#' @references \url{https://github.com/gkanogiannis/metabinR}
#'
composition_based_binning <- function(..., kMerSizeCB = 4,
                                    numOfClustersCB = 5, outputCB="CB.cluster",
                                    keepQuality = FALSE, dryRun = FALSE,
                                    gzip = FALSE, numOfThreads = 1) {
    ins <- unlist(list(...))

    composition_based_binning_checkParams(ins = ins,
                                        kMerSizeCB = kMerSizeCB,
                                        numOfClustersCB = numOfClustersCB,
                                        outputCB = outputCB,
                                        keepQuality = keepQuality,
                                        dryRun = dryRun, gzip = gzip,
                                        numOfThreads = numOfThreads)

    metatarget <- rJava::.jnew(
        class="fr/cea/ig/metatarget/MTxCB",
        class.loader = .rJava.class.loader)
    cmd <- paste("--kMerSizeCB", kMerSizeCB,
                 "--numOfClustersCB", numOfClustersCB, "--outputCB", outputCB,
                 ifelse(keepQuality, "--quality", ""),
                 ifelse(dryRun, "--dry", ""),
                 ifelse(gzip, "--gzip", ""),
                 "--numOfThreads", numOfThreads,
                 "--input", paste(ins, collapse = " --input "), sep = " ")
    ret.str <- metatarget$go(rJava::.jarray(strsplit(cmd, "\\s+")[[1]]))

    if(is.null(ret.str)) {
        return(NULL)
    } else {
        return(utils::read.table(text = ret.str, header = TRUE))
    }
}

composition_based_binning_checkParams <- function(ins, kMerSizeCB,
                                                numOfClustersCB, outputCB,
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

    if (!is.numeric(kMerSizeCB) || (is.numeric(kMerSizeCB) && kMerSizeCB<2)) {
        stop("kMerSizeCB parameter must be positive integer >1 .")
    }

    if (!is.numeric(numOfClustersCB) ||
        (is.numeric(numOfClustersCB) && numOfClustersCB<2)) {
        stop("numOfClustersCB parameter must be positive integer >1 .")
    }

    if ((!is.null(outputCB) && !methods::is(outputCB, "character")) ||
        (methods::is(outputCB, "character") && nchar(outputCB)==0)) {
        stop("outputCB must be a prefix location.")
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
