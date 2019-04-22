## The m/z grouping function 
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
#' 

mzGroup <- function(object, binGap = 0.005){
    if(class(object) != "metabCombiner")
        stop(paste(object, "is not a metabCombiner object"), sep = " ")
  
    xdata = getData(object, "x")
    ydata = getData(object, "y")
    
    if(nrow(xdata) == 0 | nrow(ydata) == 0)
        stop("missing xdata or ydata")
    
    if(!("mz" %in% names(xdata)) | !("mz" %in% names(ydata)))
        stop("m/z column missing in either xdata or ydata")
    
    object@binGap = binGap
    
    if(is.null(object@binGap) | length(object@binGap) == 0)
        stop("parameter 'binGap' must be defined")
    
    if(object@binGap <= 0 | !is.numeric(object@binGap))
        stop("parameter 'binGap' must be a positive numeric constant")
    
    if(object@binGap > 0.1)
        warning("Parameter 'binGap' is very high. Recommend value less than 0.1")
    
    xdata <- dplyr::mutate(index = 1:nrow(xdata)) %>% dplyr::arrange(mz) 
    ydata <- dplyr::mutate(index = 1:nrow(ydata)) %>% dplyr::arrange(mz) 
    
    
    ####### pick up from here!!!!!! ##########
    
    mzGrp <- rbind(dplyr::select(object@xdata@data, mz), 
                   dplyr::select(object@ydata@data, mz)) %>% 
                   dplyr::mutate(dataset = c(rep("x", nrow(object@xdata@data)),
                                rep("y", nrow(object@ydata@data)))) %>%
                   dplyr::arrange(mz)
    
    
    if(!is.numeric(mzGrp$mz) | any(is.na(mzGrp$mz)) | any(mzGrp$mz < 0)){
      stop("At least one negative, missing, or non-numeric m/z value")
    }
    
    mzGrp$groups <- .Call("binByMZ", mzGrp$mz, nrow(mzGrp), mzGrp$dataset, object@binGap)
    
    object@xdata@data$group <- dplyr::filter(mzGrp, dataset == "x")$groups 
    object@ydata@data$group <- dplyr::filter(mzGrp, dataset == "y")$groups
    
    return(object)
}
)



