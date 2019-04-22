## mzGroup 
#' @title Binning of mass spectral features in m/z dimension  
#'
#' @description The \code{mzGroup} merges m/z values of xdata &
#'     ydata into a single (sorted) list and bins the values into separate groups.
#'     Gaps between consecutive features must be less than the binGap parameter
#'     (default: 0.005) to form a group and each group must contain at least 1
#'     feature from xdata & ydata. The "group" column is then added to the "data"
#'     field for xdata & ydata
#'     
#' @param object metabCombiner object
#' 
#' @param binGap numeric. 
#' 
## 
mzGroup <- function(object, binGap = 0.005){
    if(class(object) != "metabCombiner")
        stop(base::paste(object, "is not a metabCombiner object"), sep = " ")
  
    xdata = getData(object, "x") %>% dplyr::mutate(index = 1:nrow(xdata)) 
    ydata = getData(object, "y") %>% dplyr::mutate(index = 1:nrow(ydata)) 
    
    if(nrow(xdata) == 0 | nrow(ydata) == 0)
        stop(base::paste("missing xdata or ydata in object", object, sep = " "))
    
    if(!("mz" %in% names(xdata)) | !("mz" %in% names(ydata)))
        stop("m/z column missing in either xdata or ydata")
    
    if(is.null(binGap) | length(binGap) == 0)
        stop("parameter 'binGap' must be defined")
    
    if(binGap <= 0 | !is.numeric(binGap))
        stop("parameter 'binGap' must be a positive numeric constant")
    
    if(binGap > 0.1)
        warning("parameter 'binGap' is very high. A value less than 0.1 is recommended.")
    
    object@binGap = binGap
    
    
    #forming m/z groups table
    mzGrps <- rbind(dplyr::select(xdata, mz, index), dplyr::select(ydata, mz, index)) %>% 
              dplyr::mutate(set = c(rep("x", nrow(xdata)),rep("y", nrow(ydata)))) %>%
              dplyr::arrange(mz)
    
    
    if(!is.numeric(mzGrp$mz) | any(is.na(mzGrp$mz)) | any(mzGrp$mz < 0)){
        stop("At least one negative, missing, or non-numeric m/z value")
    }
    
    mzGrps$groups <- .Call("binByMZ", 
                           mz = mzGrp$mz, 
                           length = nrow(mzGrp), 
                           datasets = mzGrp$set, 
                           gap = binGap)
    
    mzGrpsX <- dplyr::filter(mzGrp, set == "x") %>% dplyr::arrange(index)
    mzGrpsY <- dplyr::filter(mzGrp, set == "y") %>% dplyr::arrange(index)
    
    object@xdata@data$group <- mzGrpX$groups
    object@ydata@data$group <- mzGrpY$groups
    
    object = formCombinerTable(object)
    
    return(object)
}




