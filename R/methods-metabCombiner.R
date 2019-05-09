##Methods for metabCombiner object

## combinerTable
#'
#' @param object metabCombiner object.
#'
#' @returns   constructed metabCombiner table. 
## 
setMethod("combinerTable", signature = "metabCombiner", function(object){
    return(object@combinerTable)
})

## getData
#'
#' @description  
#'
#' @param object  metabCombiner object
#' 
#' @param data    Either one of 'x' or 'y'. 
#'
#' @returns If data is "x", returns data from xdata; if "y", returns data from ydata.
#'
##
setMethod("getData", signature = "metabCombiner", function(object, data = c("x", "y")){
  
    data = match.arg(data)
  
    if(data == "x")
        return(object@xdata@data)
  
    if(data == "y")
        return(object@ydata@data)
})

## getData
#'
#' @description  
#'
#' @param object  metabCombiner object
#' 
#' @param data    Either one of 'x' or 'y'. 
#'
#' @returns If data is "x", returns data from xdata; if "y", returns data from ydata.
#'
##
setMethod("getModel", signature = "metabCombiner", function(object, fit = c("loess", "gam")){
  
    fit = match.arg(fit)
  
    if(fit == "loess")
        return(object@model$loess)
  
    if(fit == "gam")
        return(object@model$gam)
})




## getSamples
#'
#' @param object metabCombiner object
#'
#' @param data Either one of 'x' or 'y'.
#'
#' @returns If data is "x", returns sample names from xdata; if "y", returns sample names from ydata.
## 
setMethod("getSamples", signature = "metabCombiner", function(object, data = c("x", "y")){
  
    data = match.arg(data)
  
    if(data == "x")
        return(object@xdata@samples)
  
    else if(data == "y")
        return(object@ydata@samples)
})

## getAnchors
#'
#' @param object metabCombiner object
#' 
#' @returns Anchor set constructed using selectAnchors().
#'  
##
setMethod("getAnchors", signature = "metabCombiner", function(object){
    return(object@anchors)
})



## plot
#'
#' @param 
#'
#' 
##
setMethod("plot", signature = "metabCombiner", function(object, fit = c("loess", "gam"),...){
    model = getModel(object, fit = fit)
    
    if(is.null(model))
        stop("missing model in metabCombiner object")
    
    table = data.frame(rtx = model$x, rty = model$y, pred = model$fitted)
    table = dplyr::arrange(table, rtx)
    
    graphics::plot(table$rtx, table$rty, type = "p", xlab = "rtx", ylab = "rty",...)
    lines(x = table$rtx, y = table$pred, col = "red", lwd = 4)
})





