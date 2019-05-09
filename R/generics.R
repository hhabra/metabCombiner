############  Generics sorted in alphabetical order #############

setGeneric("combinerTable", function(object) standardGeneric("combinerTable"))

setGeneric("getAnchors", function(object) standardGeneric("getAnchors"))

setGeneric("getModel", function(object, fit = c("loess", "gam")) standardGeneric("getModel"))

setGeneric("getData", function(object, data = c("x", "y")) standardGeneric("getData"))

setGeneric("getSamples", function(object, data = c("x", "y")) standardGeneric("getSamples"))

setGeneric("plot", function(object, ...) standardGeneric("plot"))



