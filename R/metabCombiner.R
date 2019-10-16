##
#' @title Form a metabCombiner object.
#'
#' @description
#' \code{metabCombiner()} takes two metabolomics featurelists, contained in 
#' \code{metabData} objects,and constructs a merged dataset containing groups of 
#' features detected in both featurelists with similar m/z values. 
#'
#' @param xdata metabData object. One of two datasets to be Combined.
#' 
#' @param ydata metabData object. One of two datasets to be Combined.  
#'
#' @param binGap numeric. Parameter used for grouping features by m/z. 
#' See ?mzGroup for more details.
#' 
#' @details
#' 
#' The \code{binGap} argument defines consecutive differences between
#' grouped m/z features. \code{metabCombiner} objects are used for all subsequent
#' processing steps.
#' 
#' 
#' @return a metabCombiner object 
#' 
##
metabCombiner <- function(xdata, ydata, binGap = 0.005){
    codeX = isMetabData(xdata)
    codeY = isMetabData(ydata)
    
    if(codeX)
        stop(paste("invalid xdata parameter:",combinerError(codeX, "metabData")))
    
    if(codeY)
        stop(paste("invalid ydata parameter:",combinerError(codeY, "metabData")))
    
    object <- new("metabCombiner")
    object@xdata = xdata
    object@ydata = ydata
    
    object = mzGroup(object = object, binGap = binGap)
    
    xdata = getData(object, "x") %>% dplyr::filter(group > 0) %>% 
            dplyr::arrange(group)
    
    ydata = getData(object, "y") %>% dplyr::filter(group > 0) %>% 
            dplyr::arrange(group)
    
    nGroups = max(xdata[["group"]])
    
    object@stats[["nGroups"]] = nGroups
    object = formCombinerTable(object, xdata, ydata, nGroups)
  
    return(object)
} 

