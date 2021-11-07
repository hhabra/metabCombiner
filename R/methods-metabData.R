##Methods For metabData Objects
#' @include generics.R classes.R

#' @rdname getData
#'
#' @export
##
setMethod("getData", signature = "metabData", function(object){
    data <- object@data
    return(data)
})

#' @rdname getExtra
#'
#' @export
##
setMethod("getExtra", signature = "metabData", function(object){
    extra <- object@extra
    return(extra)
})

#' @rdname getSamples
#'
#' @export
setMethod("getSamples", signature = "metabData", function(object){
    samples <- object@samples

    return(samples)
})

#' @rdname getStats
#'
#' @export
##
setMethod("getStats", signature = "metabData", function(object){
    extra <- object@stats
    return(extra)
})



##show method
setMethod("show", signature = "metabData", function(object){
    data <- getData(object)
    samples <- getSamples(object)
    extra <- getExtra(object)
    stats <- getStats(object)

    unit <- ifelse(max(data[["rt"]]) < 240, "minutes", "seconds")

    cat("A metabData object\n")
    cat("-------------------------\n")
    cat("Total Samples:", length(samples), "  Total Extra:", length(extra),"\n")

    cat("Mass Range: ", min(data[["mz"]]), "-", max(data[["mz"]]), " ",
        "m/z\n", sep = "" )
    cat("Time Range: ", min(data[["rt"]]), "-", max(data[["rt"]]), " ", unit,
        "\n", sep = "")

    cat("Input Feature Count:", stats[["input_size"]], "\n")
    cat("Final Feature Count:", stats[["final_count"]], "\n")
})


setMethod("update_md", signature = "metabData",
        function(object, data, samples, extra, stats)
{
    if(!missing(data))
        object@data <- data

    if(!missing(samples))
        object@samples <- samples

    if(!missing(extra))
        object@extra <- extra

    if(!missing(stats)){
        object@stats <- stats
    }

    return(object)
})
