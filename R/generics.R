############  Generics sorted in alphabetical order #############

setGeneric("combinerTable", function(object) standardGeneric("combinerTable"))

setGeneric("getAnchors", function(object) standardGeneric("getAnchors"))

setGeneric("getModel", function(object, fit = c("gam", "loess")) standardGeneric("getModel"))

setGeneric("getData", function(object, data = c("x", "y")) standardGeneric("getData"))

setGeneric("getExtra", function(object, data = c("x", "y")) standardGeneric("getExtra"))

setGeneric("getSamples", function(object, data = c("x", "y")) standardGeneric("getSamples"))

setGeneric("getStats", function(object) standardGeneric("getStats"))

#may want to remove this
setGeneric("plot", function(object, ...) standardGeneric("plot"))



