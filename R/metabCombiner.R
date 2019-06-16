
## metabCombiner
#'
#' @param xdata metabData object. One of two datasets to be Combined.
#' 
#' @param ydata metabData object. One of two datasets to be Combined.  
#'
#' @param binGap numeric. Parameter used for grouping features by m/z. 
#' See ?mzGroup for more details.
##
metabCombiner <- function(xdata, ydata, binGap){
    object <- new("metabCombiner")
  
    object@xdata = xdata
    object@ydata = ydata
  
    object = mzGroup(object = object, binGap = binGap)
  
    object = formCombinerTable(object)
  
    return(object)
} 

