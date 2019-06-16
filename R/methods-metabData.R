##Methods for metabData object

## Get processed metabolomics data frame from metabData object
#'
#' @param object A metabData object
#' 
#' @return Processed metabolomics feature data frame.
##
setMethod("getData", signature = "metabData", function(object){
    data = object@data
  
    return(data)  
})


## Get sample names from metabData object
#'
#' @param object A metabData object
#' 
#' @return Names of samples from metabData object.
#'
##
setMethod("getSamples", signature = "metabData", function(object){
    samples = object@samples
  
    return(samples)  
})


## Get names for "extra" columns from metabData object
#'
#' @param object A metabData object
#' 
#' @return Names of additional columns in \code{data} field of metabData
#' 
##
setMethod("getExtra", signature = "metabData", function(object){
    extra = object@extra
  
    return(extra)
})



