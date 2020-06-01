##Methods for metabCombiner objects
#' @include generics.R classes.R


#' @title Obtain metabCombiner Report
#'
#' @description Returns the constructed table reporting every possible metabolomics
#' feature pair alignment, as determined from m/z grouping step.
#'
#' @param object metabCombiner object.
#'
#' @return data.frame \code{combinedTable} merged feature alignment report. The
#' columns of the report are as follows:
#'
#' \item{idx}{Identities of features from dataset X}
#' \item{idy}{Identities of features from dataset Y}
#' \item{mzx}{m/z values of features from dataset X}
#' \item{mzy}{m/z values of features from dataset Y}
#' \item{rtx}{retention time values of features from dataset X}
#' \item{rty}{retention time values of features from dataset Y}
#' \item{rtProj}{model-projected retention time values from X to Y}
#' \item{Qx}{abundance quantile values of features from dataset X}
#' \item{Qy}{abundance quantile values of features from dataset Y}
#' \item{group}{m/z feature group}
#' \item{score}{computed similarity scores of feature pairing}
#' \item{rankX}{ranking of pairing score for X dataset features}
#' \item{rankY}{ranking of pairing score for Y dataset features}
#' \item{adductX}{adduct label of features from dataset X}
#' \item{adductY}{adduct label of features from dataset Y}
#' \item{...}{Sample and extra columns from both datasets X & Y}
#'
#' @export
setMethod("combinedTable", signature = "metabCombiner", function(object){
    return(object@combinedTable)
})

#' @title Get Ordered Retention Time Pairs
#'
#' @description
#' This returns the data frame of feature alignments used to anchor a retention
#' time projection model, constructed by \code{\link{selectAnchors}}.
#'
#' @param object metabCombiner object
#'
#' @return Data frame of anchor features
#'
#' @seealso
#' \code{\link{selectAnchors}}
#'
#' @export
setMethod("getAnchors", signature = "metabCombiner", function(object){
    return(object@anchors)
})

#' @title Obtain Last-Used Score Coefficients
#'
#' @description
#' Provides the last used weight arguments from \code{calcScores()} function.
#' Returns empty list if \code{calcScores()} has not yet been called.
#'
#' @param object  metabCombiner object
#'
#' @return A list of the last used weight parameters:
#' \item{A}{Specific weight penalizing feature m/z differences}
#' \item{B}{Specific weight penalizing relative error of retention time projection}
#' \item{C}{Specific weight penalizing differences in abundance quantiles}
#'
#' @export
setMethod("getCoefficients", signature = "metabCombiner", function(object){
    return(object@coefficients)
})


#' @title Get Fitted RT Model
#'
#' @description Returns the last fited RT projected model from a metabCombiner
#' object of type "gam" or "loess".
#'
#' @param object \code{metabCombiner} object
#'
#' @param fit Choice of either "gam" or "loess" model
#'
#' @export
setMethod("getModel", signature = "metabCombiner", function(object, fit = c("gam", "loess")){
    fit = match.arg(fit)

    if(fit == "loess")
        return(object@model[["loess"]])

    if(fit == "gam")
        return(object@model[["gam"]])
})


#' @title Get Sample Names From metabCombiner Object
#'
#' \code{metabCombiner} objects consist of two formatted metabolomics feature
#' tables. This method returns the sample names from one of the two datasets.
#'
#' @param object  metabCombiner object
#'
#' @param data   Character. One of either 'x' or 'y'.
#'
#' @return If data is "x", returns sample names for dataset X; if "y", returns
#' sample names from dataset Y.
#'
#' @export
setMethod("getSamples", signature = "metabCombiner", function(object, data = c("x", "y")){
    data = match.arg(data)

    if(data == "x")
        return(object@samples[["x"]])

    else if(data == "y")
        return(object@samples[["y"]])
})


#' @describeIn getStats Method for 'metabCombiner' object
#'
#' @export
setMethod("getStats", signature = "metabCombiner", function(object){
    return(object@stats)
})


#' @describeIn nonmatched Method for 'metabCombiner' objects
#'
#' @export
setMethod("nonmatched", signature = "metabCombiner", function(object, data = c("x", "y")){
  data = match.arg(data)

  if(data == "x")
    return(object@nonmatched[["x"]])

  if(data == "y")
    return(object@nonmatched[["y"]])
})

#' @rdname plot_Combiner
#'
#' @param x  \code{metabCombiner} object
#'
#' @param y ...
#'
#'@export
setMethod("plot", signature = "metabCombiner", function(x,y,...){
    plot_Combiner(object = x,...)
})


## Show Method
setMethod("show", signature = "metabCombiner", function(object){
    cTable = combinedTable(object)
    stats = getStats(object)

    anchors = getAnchors(object)

    coefficients = getCoefficients(object)
    A = coefficients[["A"]]
    B = coefficients[["B"]]
    C = coefficients[["C"]]

    model_gam = getModel(object)
    model_loess = getModel(object, fit = "loess")

    if(is.null(stats[["nGroups"]]))
        stats[["nGroups"]] = 0

    cat("A metabCombiner object\n")
    cat("-------------------------\n")

    if(!is.null(stats[["binGap"]]))
        cat("m/z gap size:",stats[["binGap"]], "\n")

    cat("Total Groups:", stats[["nGroups"]], "\n")

    cat("Dataset X contains", stats[["input_size_X"]], "features of which",
        stats[["grouped_size_X"]], "are grouped\n")

    cat("Dataset Y contains", stats[["input_size_Y"]], "features of which",
        stats[["grouped_size_Y"]], "are grouped\n")

    cat("combinedTable currently contains", nrow(cTable), "rows\n")

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
