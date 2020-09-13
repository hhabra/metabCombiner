##Methods for metabCombiner objects
#' @include generics.R classes.R

#' @rdname combinedTable
#'
#' @export
setMethod("combinedTable", signature = "metabCombiner", function(object){
    return(object@combinedTable)
})

#' @rdname getAnchors
#'
#' @export
setMethod("getAnchors", signature = "metabCombiner", function(object){
    return(object@anchors)
})

#' @rdname getCoefficients
#'
#' @export
setMethod("getCoefficients", signature = "metabCombiner", function(object){
    return(object@coefficients)
})

#' @rdname getModel
#'
#' @export
setMethod("getModel", function(object, fit = c("gam", "loess")){
    fit = match.arg(fit)
    if(fit == "loess")
        return(object@model[["loess"]])

    if(fit == "gam")
        return(object@model[["gam"]])
}, signature = "metabCombiner")


#' @rdname getSamples
#'
#' @export
setMethod("getSamples", function(object, data = c("x", "y")){
    data = match.arg(data)

    if(data == "x")
        return(object@samples[["x"]])

    else if(data == "y")
        return(object@samples[["y"]])
}, signature = "metabCombiner")


#' @describeIn getStats Method for 'metabCombiner' object
#'
#' @export
setMethod("getStats", signature = "metabCombiner", function(object){
    return(object@stats)
})


#' @rdname nonmatched
#'
#' @export
setMethod("nonmatched", signature = "metabCombiner",
            function(object, data = c("x", "y")){
    data = match.arg(data)

    if(data == "x")
        return(object@nonmatched[["x"]])

    if(data == "y")
        return(object@nonmatched[["y"]])
})

#' @rdname plot_fit
#'
#' @param x  \code{metabCombiner} object
#'
#' @param y ...
#'
#'@export
setMethod("plot", signature = "metabCombiner", function(x,y,...){
    plot_fit(object = x,...)
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
        cat("Last Scoring Weights Used:\nm/z weight = ", A, ", rt weight = ",
            B,  ", quantile weight = ", C, sep = "")
    return(invisible())
})
