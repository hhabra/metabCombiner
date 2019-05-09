##Methods for metabData object

## getData
#'
#'
#'
##
setMethod("getData", signature = "metabData", function(object){
    data = object@data
  
    return(data)  
})


## getSamples
#'
#'
#'
##
setMethod("getSamples", signature = "metabData", function(object){
    samples = object@samples
  
    return(samples)  
})





