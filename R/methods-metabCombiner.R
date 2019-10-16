##Methods for metabCombiner object

#' @title Obtain Combiner Feature Alignment Report
#'
#' @description Obtain constructed table reporting every possible metabolomics
#' feature alignment (as determined by m/z binGap parameter in metabCombiner().
#'
#' @param object metabCombiner object.
#'
#' @return Data frame combinerTable merged feature alignment report. The columns
#' of the report are as follows:
#'
#' \item{idx}{Identities of features from dataset X}
#' \item{idy}{Identities of features from dataset Y}
#' \item{mzx}{m/z values of features from dataset X}
#' \item{mzy}{m/z values of features from dataset Y}
#' \item{rtx}{retention time values of features from dataset X}
#' \item{rty}{retention time values of features from dataset Y}
#' \item{rtProj}{model-projected (X->Y) retention times values}
#' \item{Qx}{abundance quantile values of features from dataset X}
#' \item{Qy}{abundance quantile values of features from dataset Y}
#' \item{group}{m/z feature group of feature pairing}
#' \item{score}{computed similarity scores of feature pairing}
#' \item{rankX}{ranking of pairing score for X dataset features}
#' \item{rankY}{ranking of pairing score for Y dataset features}
#' \item{adductX}{adduct label of features from dataset X}
#' \item{adductY}{adduct label of features from dataset Y}
#' \item{...}{Sample and extra columns from both datasets X & Y}
#'
#' @exportMethod
##
setMethod("combinerTable", signature = "metabCombiner", function(object){
    return(object@combinerTable)
})

#' @title Get Ordered Retention Time Pairs
#'
#' @param object metabCombiner object
#'
#' @return anchor set constructed using selectAnchors().
#'
#' @exportMethod
##
setMethod("getAnchors", signature = "metabCombiner", function(object){
    return(object@anchors)
})

#' @title Obtain last-used set of score coefficients
#'
#' @description
#' Provides the last used weight arguments from \code{calcScores()} function.
#' Returns NULL if \code{calcScores()} has not yet been called.
#'
#' @param object  metabCombiner object
#'
#' @return A list of the last used weight parameters:
#' \item{A}{Specific weight penalizing feature m/z differences}
#' \item{B}{Specific weight penalizing relative error of retention time projection}
#' \item{C}{Specific weight penalizing differences in abundance quantiles}
#'
#' @exportMethod
##
setMethod("getCoefficients", signature = "metabCombiner", function(object){
    return(object@coefficients)
})


#' @title Get Datasets From metabCombiner Object.
#'
#' @description
#'
#' @param object  metabCombiner object
#'
#' @param data    Either one of 'x' or 'y'.
#'
#' @return If data is "x", returns x dataset ; if "y", returns y dataset.
#'
#' @exportMethod
##
setMethod("getData", signature = "metabCombiner", function(object, data = c("x", "y")){
    data = match.arg(data)

    if(data == "x")
        return(object@xdata@data)

    if(data == "y")
        return(object@ydata@data)
})

#' @title Get Fitted RT Model
#'
#' @description
#'
#' @param object  metabCombiner object
#'
#' @param fit   Choice of model. Either "loess" or "gam".
#'
#' @return Nonlinear retention time fitting model.
#'
#' @exportMethod
##
setMethod("getModel", signature = "metabCombiner", function(object, fit = c("gam", "loess")){
    fit = match.arg(fit)

    if(fit == "loess")
        return(object@model[["loess"]])

    if(fit == "gam")
        return(object@model[["gam"]])
})


#' @title Get Sample Names From metabCombiner Object
#'
#' @param object  metabCombiner object
#'
#' @param data    Either one of 'x' or 'y'.
#'
#' @return If data is "x", returns sample names from xdata; if "y", returns
#' sample names from ydata.
#'
#' @exportMethod
##
setMethod("getSamples", signature = "metabCombiner", function(object, data = c("x", "y")){

    data = match.arg(data)

    if(data == "x")
        return(object@xdata@samples)

    else if(data == "y")
        return(object@ydata@samples)
})


#' @title Get metabCombiner Statistics
#'
#' @param object  metabCombiner object
#'
#' @return A list containing statistics from metabCombiner object.
#'
#' @exportMethod
##
setMethod("getStats", signature = "metabCombiner", function(object){
    return(object@stats)
})


#' @title Plot nonlinear fit between retention times
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
        stop(paste("object missing model of type", fit))

    if(fit == "loess")
        data = data.frame(rtx = model[["x"]],
                          rty = model[["y"]],
                          pred = model[["fitted"]],
                          weights = model[["weights"]])

    else if (fit == "gam")
        data = data.frame(rtx = model$model$rtx,
                          rty = model$model$rty,
                          pred = model$fitted.values,
                          weights = model$prior.weights)

    data = dplyr::arrange(data, rtx)

    if(missing(pcol))
        pcol = "black"
    if(missing(lcol))
        lcol = "red"
    if(missing(lwd))
        lwd = 4

    graphics::plot(data[["rtx"]], data[["rty"]], type = "p", col = pcol,...)
    graphics::lines(x = data[["rtx"]], y = data[["pred"]],
                    col = lcol, lwd = lwd,...)
})


## Show Method
setMethod("show", signature = "metabCombiner", function(object){
    xdata = getData(object)
    ydata = getData(object, data = "y")
    cTable = combinerTable(object)

    stats = getStats(object)

    anchors = getAnchors(object)

    coefficients = getCoefficients(object)
    A = coefficients[["A"]]
    B = coefficients[["B"]]
    C = coefficients[["C"]]

    model_gam = getModel(object)
    model_loess = getModel(object, fit = "loess")

    groupedX = 0
    groupedY = 0

    if(nrow(xdata) > 0)
        groupedX = sum(xdata[["group"]] != 0)

    if(nrow(ydata) > 0)
        groupedY = sum(ydata[["group"]] != 0)

    if(is.null(stats[["nGroups"]]))
        stats[["nGroups"]] = 0

    cat("A metabCombiner object\n")
    cat("-------------------------\n")
    cat("Gap Parameter:", object@binGap, "\n")

    cat("Total Groups:", stats[["nGroups"]], "\n")
    if(!is.null(stats[["binGap"]]))
        cat("m/z gap size:",stats[["binGap"]], "\n")

    cat("Dataset X contains", nrow(xdata), "features of which", groupedX,
        "are grouped\n")

    cat("Dataset Y contains", nrow(ydata), "features of which", groupedY,
        "are grouped\n")

    cat("combinerTable currently contains", nrow(cTable), "rows\n")

    if(nrow(anchors) > 0)
        cat(nrow(anchors), "ordered pairs selected to fit RT model\n")

    if(!is.null(model_gam) & !is.null(stats[["best_k"]]))
        cat("Object contains a GAM model with k =", stats[["best_k"]], "\n")

    if(!is.null(model_loess) & !is.null(stats[["best_span"]]))
        cat("Object contains a LOESS model with span =", stats[["best_span"]],
            "\n")

    if(!is.null(A) & !is.null(B) & !is.null(C))
        cat("Last Scoring Parameters Used:\nm/z weight = ", A, ", rt weight = ",
            B,  ", quantile weight = ", C, sep = "")

    return(invisible())
})


