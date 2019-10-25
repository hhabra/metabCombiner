#' @title Print a Combiner Report to file.
#'
#' @description Prints a \code{combinerTable} report to a file, specified by
#'              \code{file} argument. The output will contain an empty line
#'              separating each m/z group.
#'
#' @param object  metabCombiner object or combinerTable.
#'
#' @param file  character string naming an output file path.
#'
#' @param sep   Character field separator. Values within each row are separated
#'              by this string.
#'
#' @export
##
write2file <- function(object, file = "", sep = ","){

    if(isCombinerTable(object) == 0)
        cTable = object

    else{
        code = isMetabCombiner(object)

        if(code)
            stop(combinerError(code, "metabCombiner"))

        cTable = combinerTable(object)
    }

    sep = as.character(sep)
    if(base::nchar(sep) > 1)
        stop("parameter 'sep' must be a length <= 1 character")

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
                        PACKAGE = "Combiner"
          )

    return(invisible())
}


