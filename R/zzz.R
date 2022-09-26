#' @importFrom utils citation
.onLoad <- function(libname, pkgname) {

    rJava::.jpackage(
        name = pkgname,
        lib.loc=libname,
        own.loader = TRUE,
        parameters = c(getOption("java.parameters")
                        ,"-Djava.awt.headless=true"
                        ,"-XX:+UseG1GC"
                        ,"-XX:+UseStringDeduplication"
                        )
    )
}

.onAttach <- function(libname, pkgname){
    requireNamespace("utils")
    cit<-citation(pkgname)
    txt<-paste(c(format(cit,"citation")),collapse="\n\n")
    packageStartupMessage(txt)
}
