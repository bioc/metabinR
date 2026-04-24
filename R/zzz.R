#' JVM options used when metabinR loads
#'
#' Returns the vector of JVM flags passed to \code{\link[rJava]{.jpackage}}
#' on package load. Set \code{options(metabinR.jvm.flags = c(...))} before
#' loading the package to override; set \code{options(java.parameters = ...)}
#' to prepend heap-size flags (e.g. \code{"-Xmx4g"}) in the usual rJava way.
#'
#' @return A character vector of JVM flags.
#' @examples
#' metabinR_jvm_options()
#' @export
metabinR_jvm_options <- function() {
    user_flags <- getOption("metabinR.jvm.flags", default = NULL)
    defaults <- c(
        "-Djava.awt.headless=true",
        "-XX:+UseG1GC",
        "-XX:+UseStringDeduplication"
    )
    c(getOption("java.parameters"),
      if (is.null(user_flags)) defaults else user_flags)
}

.onLoad <- function(libname, pkgname) {
    rJava::.jpackage(
        name = pkgname,
        lib.loc = libname,
        own.loader = TRUE,
        parameters = metabinR_jvm_options()
    )
    # devtools::load_all() does not always wire the own.loader classpath;
    # add every JAR under inst/java / java explicitly as a fallback.
    jar_dir <- system.file("java", package = pkgname, lib.loc = libname)
    jars <- if (nzchar(jar_dir)) {
        list.files(jar_dir, pattern = "\\.jar$", full.names = TRUE)
    } else character(0)
    loader <- tryCatch(get(".rJava.class.loader", envir = asNamespace(pkgname)),
                       error = function(e) NULL)
    for (jar in jars) {
        if (!is.null(loader)) {
            try(rJava::.jcall(loader, "V", "addClassPath", jar), silent = TRUE)
        }
        try(rJava::.jaddClassPath(jar), silent = TRUE)
    }
}
