##Methods for metabCombiner object

## Get combinerTable from metabCombiner object.
#'
#' @param object metabCombiner object.
#'
#' @return   constructed metabCombiner table. 
## 
setMethod("combinerTable", signature = "metabCombiner", function(object){
    return(object@combinerTable)
})

## Get 'x' or 'y' dataset from metabCombiner object.
#'
#' @description  
#'
#' @param object  metabCombiner object
#' 
#' @param data    Either one of 'x' or 'y'. 
#'
#' @return If data is "x", returns x dataset ; if "y", returns y dataset.
#'
##
setMethod("getData", signature = "metabCombiner", function(object, data = c("x", "y")){
    data = match.arg(data)
  
    if(data == "x")
        return(object@xdata@data)
  
    if(data == "y")
        return(object@ydata@data)
})

## Get fitted regression model from metabCombiner object.
#'
#' @description  
#'
#' @param object  metabCombiner object
#' 
#' @param fit   Choice of model. Either "loess" or "gam". 
#'
#' @return Nonlinear retention time fitting model.
#'
##
setMethod("getModel", signature = "metabCombiner", function(object, fit = c("loess", "gam")){
  
    fit = match.arg(fit)
  
    if(fit == "loess")
        return(object@model$loess)
  
    if(fit == "gam")
        return(object@model$gam)
})


## Get 'x' or 'y' sample names from metabCombiner object.
#'
#' @param object  metabCombiner object
#' 
#' @param data    Either one of 'x' or 'y'.
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

## Get list of "anchor" ordered retention time pairs
#'
#' @param object metabCombiner object
#' 
#' @returns anchor set constructed using selectAnchors().
#'  
##
setMethod("getAnchors", signature = "metabCombiner", function(object){
    return(object@anchors)
})



## plot nonlinear retention time curve 
#'
#' @param object metabCombiner object
#'
#' @param fit Choice of model. Either "loess" or "gam".
#' 
#' @param ... Other variables passed into graphics::plot
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