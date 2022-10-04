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
#' @param ... Input fasta/fastq files locations
#'     (uncompressed or gzip compressed).
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
#' @param numOfThreads Number of threads to use.
#'
#' @return A \code{\link[base]{data.frame}} of the binning assignments.
#'     Return value contains \code{numOfClustersAB + 2} columns.
#' \itemize{
#'     \item \code{read_id} : read identifier from fasta header
#'     \item \code{ABxCB} : read was assigned to this ABxCB cluster index
#'     \item \code{ABxCB.n} : read to cluster ABxCB.n distance
#' }
#' @export
#'
#' @examples
#' hierarchical_binning(
#'     system.file("extdata", "reads.metagenome.fasta.gz",package = "metabinR"),
#'     dryRun = TRUE, kMerSizeAB = 4, kMerSizeCB = 2
#' )
#' @author Anestis Gkanogiannis, \email{anestis@@gkanogiannis.com}
#' @references \url{https://github.com/gkanogiannis/metabinR}
#'

hierarchical_binning <- function(..., eMin = 1, eMax = 0, kMerSizeAB = 10,
                                    kMerSizeCB = 4, genomeSize = 3000000,
                                    numOfClustersAB =3, outputC="ABxCB.cluster",
                                    keepQuality = FALSE, dryRun = FALSE,
                                    gzip = FALSE, numOfThreads = 1) {
    ins <- unlist(list(...))

    hierarchical_binning_checkParams(ins = ins, eMin = eMin, eMax = eMax,
                                            kMerSizeAB = kMerSizeAB,
                                            kMerSizeCB = kMerSizeCB,
                                            genomeSize = genomeSize,
                                            numOfClustersAB = numOfClustersAB,
                                            outputC = outputC,
                                            keepQuality = keepQuality,
                                            dryRun = dryRun, gzip = gzip,
                                            numOfThreads = numOfThreads)

    metatarget <- rJava::.jnew(
        class="fr/cea/ig/metatarget/MTxABxCB",
        class.loader = .rJava.class.loader)
    cmd <- paste("--eMin", eMin, "--eMax", eMax, "--kMerSizeAB", kMerSizeAB,
                 "--kMerSizeCB", kMerSizeCB, "--genomeSize", genomeSize,
                 "--numOfClustersAB", numOfClustersAB, "--outputC", outputC,
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

hierarchical_binning_checkParams <- function(ins, eMin, eMax, kMerSizeAB,
                                                    kMerSizeCB, genomeSize,
                                                    numOfClustersAB, outputC,
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

    if (!is.numeric(kMerSizeCB) || (is.numeric(kMerSizeCB) && kMerSizeCB<2)) {
        stop("kMerSizeCB parameter must be positive integer >1 .")
    }

    if (!is.numeric(genomeSize) ||
        (is.numeric(genomeSize) && genomeSize<1)) {
        stop("genomeSize parameter must be positive integer.")
    }

    if (!is.numeric(numOfClustersAB) ||
        (is.numeric(numOfClustersAB) && numOfClustersAB<2)) {
        stop("numOfClustersAB parameter must be positive integer >1 .")
    }

    if ((!is.null(outputC) && !methods::is(outputC, "character")) ||
        (methods::is(outputC, "character") && nchar(outputC)==0)) {
        stop("outputC must be a prefix location.")
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
