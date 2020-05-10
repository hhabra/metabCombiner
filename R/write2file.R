#' @title Print metabCombiner Report to File.
#'
#' @description Prints a \code{combinedTable} report to a file, specified by
#'              \code{file} argument. The output will have an empty line between
#'              each separate m/z group for ease of viewing.
#'
#' @param object  \code{metabCombiner} object or \code{combinedTable}.
#'
#' @param file  character string naming the output file path
#'
#' @param sep   Character field separator. Values within each row are separated
#'              by this character.
#'
#' @examples
#' \dontrun{
#' library(metabCombiner)
#' data(plasma30)
#' data(plasma20)
#'
#' p30 <- metabData(plasma30, samples = "CHEAR")
#' p20 <- metabData(plasma20, samples = "Red", rtmax = 17.25)
#' p.combined = metabCombiner(xdata = p30, ydata = p20, binGap = 0.0075)
#'
#' p.combined = selectAnchors(p.combined, tolMZ = 0.003, tolQ = 0.3, windY = 0.02)
#' p.combined = fit_gam(p.combined, k = seq(12,20,2), iterFilter = 1)
#' p.combined = calcScores(p.combined, A = 90, B = 14, C = 0.5)
#'
#' ###using metabCombiner object as input
#' write2file(p.combined, file = "plasma-combined.csv", sep = ",")
#'
#' ###using combinedTable report as input
#' cTable = combinedTable(p.combined)
#' write2file(cTable, file = "plasma-combined.txt", sep = "\t")
#' }
#' @export
##
write2file <- function(object, file = "", sep = ","){
    if(iscombinedTable(object) == 0)
        cTable = object

    else{
        code = isMetabCombiner(object)

        if(code)
            stop(combinerError(code, "metabCombiner"))

        cTable = combinedTable(object)
    }

    sep = as.character(sep)
    if(base::nchar(sep) > 1)
        stop("argument 'sep' must be a length <= 1 character")

    if(sep == ".")
        stop("invalid input for argument 'sep'")

    #bugfix for character columns containing sep character
    colclasses = as.character(sapply(cTable, class))
    ccs = which(colclasses == "character")
    cTable[,ccs] <- lapply(cTable[,ccs], function(x) gsub(sep,".",x))

    if(!S4Vectors::isSorted(cTable[["group"]])){
        cTable = cTable[with(cTable, order(group, desc(score))), ]
    }

    if(any(is.na(cTable[["group"]])))
        stop("Missing value detected in group column")

    group = cTable[["group"]]

    ##print the first line, plus column names
    write.table(cTable[1,],
                file = file,
                row.names = FALSE,
                sep = sep,
                na = "")

    lines = apply(cTable,1, function(row) base::paste(row, collapse = sep))

    .Call("write2file", lines = lines,
                        file = file,
                        groups = group,
                        PACKAGE = "metabCombiner"
          )

    return(invisible())
}


