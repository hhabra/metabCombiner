
.onLoad <- function(libname, pkgname) {
    quietly <- getOption('quietly')
    options(quietly = T)
    suppressPackageStartupMessages(
        requireNamespace("dplyr", quietly = TRUE)
    )
}


.onUnload <- function(libpath) {
    library.dynam.unload("metabCombiner", libpath)
}


