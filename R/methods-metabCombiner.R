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
setMethod("getModel", signature = "metabCombiner", function(object, fit = c("gam", "loess")){
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
#' @return If data is "x", returns sample names from xdata; if "y", returns sample names from ydata.
## 
setMethod("getSamples", signature = "metabCombiner", function(object, data = c("x", "y")){
  
    data = match.arg(data)
  
    if(data == "x")
        return(object@xdata@samples)
  
    else if(data == "y")
        return(object@ydata@samples)
})


## Get statistics field of metabCombiner object.
#'
#' @param object  metabCombiner object
#' 
#' @return A list containing statistics from metabCombiner object.
##
setMethod("getStats", signature = "metabCombiner", function(object){
    return(object@stats)
})


## plot nonlinear retention time fitting curve between feature tables
#'
#' @param object metabCombiner object
#'
#' @param fit Choice of model. Either "gam" or "loess".
#' 
#' @param pcol color of the points (ordered pairs) in the plot.
#' 
#' @param lcol color of the fitted line in the plot
#' 
#' @param ... Other variables passed into graphics::plot
#' 
##
setMethod("plot", signature = "metabCombiner", 
          function(object, fit = c("gam","loess"), pcol, lcol, lwd,...){
    
    fit = match.arg(fit)  
    model = getModel(object, fit = fit)
    
    if(is.null(model))
        stop("missing model in metabCombiner object")
    
    if(fit == "loess")
        data = data.frame(rtx = model$x, rty = model$y, pred = model$fitted,
                           weights = model$weights)
    else if (fit == "gam")
        data = data.frame(rtx = model$model$rtx, rty = model$model$rty,
                           pred = model$fitted.values, weights = model$prior.weights)
    
    data = dplyr::arrange(data, rtx)
    
    if(missing(pcol))
        pcol = "black"
    if(missing(lcol))
        lcol = "red"
    if(missing(lwd))
        lwd = 4
      
    graphics::plot(data$rtx, data$rty, type = "p", col = pcol,...)
    lines(x = data$rtx, y = data$pred, col = lcol, lwd = lwd,...)
})