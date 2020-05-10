##Methods For metabData Objects
#' @include generics.R classes.R


#' @describeIn getData Method for "metabData" objects
#'
#' @export
##
setMethod("getData", signature = "metabData", function(object){
    data = object@data
    return(data)
})


#' @title Obtain Metabolomics Sample Names
#'
#' @param object A metabData object
#'
#' @return Names of samples of formatted dataset.
#'
#' @export
setMethod("getSamples", signature = "metabData", function(object){
    samples = object@samples

    return(samples)
})


#' @describeIn getExtra Method for "metabData" Objects
#'
#' @export
##
setMethod("getExtra", signature = "metabData", function(object){
    extra = object@extra

    return(extra)
})


##
#' @describeIn getStats Method for 'metabData' object
#'
#' @export
##
setMethod("getStats", signature = "metabData", function(object){
    extra = object@stats

    return(extra)
})



##show method
setMethod("show", signature = "metabData", function(object){
    data = getData(object)
    samples = getSamples(object)
    stats = getStats(object)

    timeUnit = ifelse(max(data[["rt"]]) < 180, "minutes", "seconds")

    cat("A metabData object\n")
    cat("-------------------------\n")
    cat("Total Samples:", length(samples), "\n")

    cat("Mass Range: ", min(data[["mz"]]), "-", max(data[["mz"]]), " ",
        "m/z\n", sep = "" )
    cat("Time Range: ", min(data[["rt"]]), "-", max(data[["rt"]]), " ", timeUnit,
        "\n", sep = "")

    cat("Input Feature Count:", stats[["input_size"]], "\n")
    cat("Final Feature Count:", stats[["final_count"]], "\n")
})


