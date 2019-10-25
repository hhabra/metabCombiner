##Methods For metabData Objects
#' @include generics.R classes.R


#' @title Obtain processed metabolomics data frame
#'
#' @description
#' Extracts the formatted feature table contained within metabData object.
#'
#' @param object A metabData object
#'
#' @return Processed metabolomics feature data frame.
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
##
setMethod("getSamples", signature = "metabData", function(object){
    samples = object@samples

    return(samples)
})


#' @title Get Names of Additional Data Fields
#'
#' @param object A metabData object
#'
#' @return Names of 'extra' columns in \code{data} field of metabData
#'
#' @export
##
setMethod("getExtra", signature = "metabData", function(object){
    extra = object@extra

    return(extra)
})


##
#' @title Obtain metabData Statistics
#'
#' @param object A metabData object
#'
#' @return Feature statistics from metabData object
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


