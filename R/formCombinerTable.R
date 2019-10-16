## 
#' @title Form Combiner Report Table
#'  
#' @description  Takes previously computed m/z groups using \code{mzGroup()}
#' and creates merged \code{combinerTable} consisting of all possible feature
#' alignments, all initialized with equivalent score and ranking.    
#'  
#' @param object metabCombiner object
#' 
#' @param xdata data frame. A processed metabolomics feature table.
#' 
#' @param ydata data frame. A processed metabolomics feature table.
#' 
#' @param nGroups integer. Total number of computed feature groups
#' 
#' @return metabCombiner object with initialized combinerTable data frame.
## 
formCombinerTable <- function(object, xdata, ydata, nGroups){
    if(class(object) != "metabCombiner")
        stop(paste(object, "is not a metabCombiner object"))
  
    groupCountX = as.integer(table(xdata[["group"]]))
    groupCountY = as.integer(table(ydata[["group"]]))
    
    if(any(groupCountX * groupCountY > 10000))
        stop("Irregular group size detected (n > 10000)! Check m/z values.")
    
    totalRows = sum(groupCountX * groupCountY)
    
    xreps = rep(groupCountY, times = groupCountX)
    xCombine = dplyr::slice(xdata, rep(1:n(), times = xreps))
    
    yreps = lapply(1:nGroups, function(number){
        counts = base::rep(x = which(ydata[["group"]] == number), 
                       times = groupCountX[number])
        return(counts)
    }) %>% unlist()
    
    yCombine = ydata[yreps,]
    
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
                        xCombine[,7:ncol(xCombine)],
                        yCombine[,7:ncol(yCombine)],
                        stringsAsFactors = FALSE, check.names = FALSE
    )
                        
  
    object@combinerTable = cTable
  
    return(object)
}












