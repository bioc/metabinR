#' @importFrom rJava .jnew .jcall .jarray .jlong
#' @importFrom S4Vectors DataFrame
#' @importFrom utils read.table
#' @importFrom cli cli_abort
NULL

# JNI return-type stub used by the typed run(...) methods.
.metabinR_JNI_RET <- "Ljava/lang/String;"

# JNI signatures for the typed run(...) methods on each Java entry class.
.metabinR_sig_AB <- paste0(
    "([Ljava/lang/String;IIIILjava/lang/String;ZZZI)",
    .metabinR_JNI_RET
)
.metabinR_sig_CB <- paste0(
    "([Ljava/lang/String;IILjava/lang/String;ZZZI)",
    .metabinR_JNI_RET
)
.metabinR_sig_ABxCB <- paste0(
    "([Ljava/lang/String;IIIIIIJLjava/lang/String;ZZZI)",
    .metabinR_JNI_RET
)

.metabinR_jclass <- c(
    AB    = "fr/cea/ig/metatarget/MTxAB",
    CB    = "fr/cea/ig/metatarget/MTxCB",
    ABxCB = "fr/cea/ig/metatarget/MTxABxCB"
)

# Build a MetabinResult from tab-separated text produced by the Java run().
.parse_assignments <- function(ret_str, inputs, params, algo) {
    if (is.null(ret_str) || !nzchar(ret_str)) {
        cli::cli_abort(
            "Java backend returned no assignments (algorithm: {algo}).",
            class = "metabinR_error_empty_result"
        )
    }
    df <- utils::read.table(
        text = ret_str, header = TRUE, sep = "\t",
        stringsAsFactors = FALSE, check.names = FALSE
    )
    new("MetabinResult",
        assignments = S4Vectors::DataFrame(df, check.names = FALSE),
        parameters  = params,
        inputs      = inputs,
        algorithm   = algo)
}

# Dispatch to the right typed Java entry point.
.call_bridge <- function(algo, inputs, params) {
    jclass <- .metabinR_jclass[[algo]]
    loader <- tryCatch(
        get(".rJava.class.loader", envir = asNamespace("metabinR")),
        error = function(e) NULL
    )
    obj <- if (is.null(loader)) {
        rJava::.jnew(class = jclass)
    } else {
        rJava::.jnew(class = jclass, class.loader = loader)
    }
    jinputs <- rJava::.jarray(as.character(inputs))

    ret_str <- switch(
        algo,
        AB = rJava::.jcall(
            obj, .metabinR_JNI_RET, "run",
            jinputs,
            as.integer(params$eMin), as.integer(params$eMax),
            as.integer(params$kMerSizeAB), as.integer(params$numOfClustersAB),
            as.character(params$outputAB),
            as.logical(params$keepQuality), as.logical(params$dryRun),
            as.logical(params$gzip), as.integer(params$numOfThreads)
        ),
        CB = rJava::.jcall(
            obj, .metabinR_JNI_RET, "run",
            jinputs,
            as.integer(params$kMerSizeCB),
            as.integer(params$numOfClustersCB),
            as.character(params$outputCB),
            as.logical(params$keepQuality), as.logical(params$dryRun),
            as.logical(params$gzip), as.integer(params$numOfThreads)
        ),
        ABxCB = rJava::.jcall(
            obj, .metabinR_JNI_RET, "run",
            jinputs,
            as.integer(params$eMin), as.integer(params$eMax),
            as.integer(params$kMerSizeAB), as.integer(params$kMerSizeCB),
            as.integer(params$numOfClustersAB),
            as.integer(params$readLength %||% 0L),
            rJava::.jlong(as.numeric(params$genomeSize)),
            as.character(params$outputC),
            as.logical(params$keepQuality), as.logical(params$dryRun),
            as.logical(params$gzip), as.integer(params$numOfThreads)
        ),
        cli::cli_abort(
            "Unknown algorithm {.val {algo}}.",
            class = "metabinR_error_bad_algo"
        )
    )

    .parse_assignments(ret_str, inputs = inputs, params = params, algo = algo)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Accept a list of inputs from `...` and resolve each element to a readable
# file path. Character scalars/vectors pass through; in-memory sequence
# objects (DNAStringSet / QualityScaledDNAStringSet / ShortReadQ) are staged
# to a tempfile so the Java backend — which only reads from disk — can
# consume them.
.resolve_inputs <- function(args) {
    if (length(args) == 0L) {
        cli::cli_abort(
            "No input sequences or FASTA/FASTQ files were provided.",
            class = "metabinR_error_no_input"
        )
    }
    out <- character(0)
    for (a in args) {
        if (is.null(a)) next
        if (is.character(a)) {
            out <- c(out, a)
        } else if (methods::is(a, "ShortReadQ")) {
            path <- tempfile(fileext = ".fastq")
            ShortRead::writeFastq(a, file = path, compress = FALSE)
            out <- c(out, path)
        } else if (methods::is(a, "QualityScaledDNAStringSet") ||
                   methods::is(a, "DNAStringSet")) {
            path <- tempfile(fileext = ".fasta")
            Biostrings::writeXStringSet(a, filepath = path, format = "fasta",
                                        compress = FALSE)
            out <- c(out, path)
        } else {
            cli::cli_abort(
                c("Unsupported input type: {.cls {class(a)[1]}}.",
                  "i" = paste0("Pass character file paths, a DNAStringSet, ",
                               "or a ShortReadQ object.")),
                class = "metabinR_error_input_type"
            )
        }
    }
    if (length(out) == 0L) {
        cli::cli_abort(
            "No input sequences or FASTA/FASTQ files were provided.",
            class = "metabinR_error_no_input"
        )
    }
    out
}
