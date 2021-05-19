#' @title Print metabCombiner Report to File.
#'
#' @description Prints a \code{combinedTable} report to a file, specified by
#'              \code{file} argument. Output file has an empty line between
#'              each separate m/z group for ease of viewing.
#'
#' @param object  \code{metabCombiner} object or \code{combinedTable}
#'
#' @param file character string naming the output file path
#'
#' @param sep   Character field separator. Values within each row are separated
#'              by this character.
#'
#' @return no values returned
#'
#' @examples
#' data(plasma30)
#' data(plasma20)
#'
#' p30 <- metabData(plasma30, samples = "CHEAR")
#' p20 <- metabData(plasma20, samples = "Red", rtmax = 17.25)
#' p.comb <- metabCombiner(xdata = p30, ydata = p20, binGap = 0.0075)
#'
#' p.comb <- selectAnchors(p.comb, tolmz = 0.003, tolrtq = 0.3, windy = 0.02)
#' p.comb <- fit_gam(p.comb, k = 20, iterFilter = 1)
#' p.comb <- calcScores(p.comb, A = 90, B = 14, C = 0.5)
#' p.comb <- labelRows(p.comb, maxRankX = 2, maxRankY = 2, remove = TRUE)
#'
#' \donttest{
#' ###using metabCombiner object as input
#' write2file(p.comb, file = "plasma-combined.csv", sep = ",")
#'
#' ###using combinedTable report and feature data as input
#' cTable <- cbind.data.frame(combinedTable(p.comb), featdata(p.comb))
#' write2file(cTable, file = "plasma-combined.txt", sep = "\t")
#' }
#' @export
##
write2file <- function(object, file, sep = ","){
    if(isCombinedTable(object) == 0) cTable = object
    else if(isMetabCombiner(object) == 0)
        cTable = cbind.data.frame(combinedTable(object), featdata(object))
    else{
        combinerCheck(isMetabCombiner(object), "metabCombiner", "warning")
        combinerCheck(isCombinedTable(object), "combinedTable", "warning")
        stop("input object is not a valid metabCombiner or combinedTable")
    }
    sep <- as.character(sep)
    if(base::nchar(sep) > 1)
        stop("argument 'sep' must be a length <= 1 character")

    if(sep == "." | sep == ";")
        stop("use of ; or . forbidden for sep argument")

    #bugfix for character columns containing sep character
    colclasses <- as.character(lapply(cTable, class))
    ccs <- which(colclasses == "character")
    cTable[,ccs] <- lapply(cTable[,ccs], function(x) gsub(sep,".",x))

    if(!S4Vectors::isSorted(cTable[["group"]]))
        cTable <- cTable[with(cTable, order(group, desc(score))), ]

    if(any(is.na(cTable[["group"]])))
        stop("missing value detected in group column")

    group <- cTable[["group"]]
    write.table(cTable[1,], file = file, row.names = FALSE, sep = sep, na = "")

    if(!file.exists(file))
        stop("invalid file argument- cannot write to indicated file")

    if(any(duplicated(names(cTable))))
        lines <- do.call(paste, c(cTable, list(sep = sep)))
    else
        lines <- tidyr::unite(cTable, col = "dat", names(cTable), sep = ",")$dat

    .Call("write2file", lines = lines, file = file, groups = group, code = 0,
            PACKAGE = "metabCombiner")

    return()
}


