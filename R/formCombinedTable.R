##
#' @title Form Combiner Report Table
#'
#' @description  Takes previously computed m/z groups using \code{mzGroup()}
#' and creates merged \code{combinedTable} consisting of all possible feature
#' alignments, all initialized with equivalent score and ranking.
#'
#' @param object metabCombiner object
#'
#' @param xset data frame. A processed metabolomics feature table.
#'
#' @param yset data frame. A processed metabolomics feature table.
#'
#' @param nGroups integer. Total number of computed feature groups
#'
#' @return metabCombiner object with initialized combinedTable data frame.
##
formCombinedTable <- function(object, xset, yset, nGroups){
    if(class(object) != "metabCombiner")
        stop(paste(object, "is not a metabCombiner object"))

    groupCountX = as.integer(table(xset[["group"]]))
    groupCountY = as.integer(table(yset[["group"]]))

    if(any(groupCountX * groupCountY > 10000))
        stop("Irregular group size detected (n > 10000)! Check m/z values.")

    totalRows = sum(groupCountX * groupCountY)

    xreps = rep(groupCountY, times = groupCountX)
    xCombine = dplyr::slice(xset, rep(1:n(), times = xreps))

    yreps = lapply(1:nGroups, function(number){
        counts = base::rep(x = which(yset[["group"]] == number),
                           times = groupCountX[number])
        return(counts)
    }) %>% unlist()

    yCombine = yset[yreps,]

    #combine groups into data frame
    cTable = data.frame(idx = xCombine[["id"]], idy = yCombine[["id"]],
                        mzx = xCombine[["mz"]], mzy = yCombine[["mz"]],
                        rtx = xCombine[["rt"]], rty = yCombine[["rt"]],
                        rtProj = numeric(totalRows),
                        Qx = xCombine[["Q"]], Qy = yCombine[["Q"]],
                        group = xCombine[["group"]],
                        score = rep(1, totalRows),
                        rankX = as.integer(1),
                        rankY = as.integer(1),
                        adductx = xCombine[["adduct"]],
                        adducty = yCombine[["adduct"]],
                        xCombine[object@samples[["x"]]],
                        xCombine[object@extra[["x"]]],
                        yCombine[object@samples[["y"]]],
                        yCombine[object@extra[["y"]]],
                        stringsAsFactors = FALSE, check.names = FALSE
    )

    object@combinedTable = cTable

    return(object)
}
