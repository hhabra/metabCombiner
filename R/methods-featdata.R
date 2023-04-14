##Methods for extracting feature meta-data
#' @include generics.R classes.R


#' @rdname adductdata
#'
#' @export
setMethod("adductData",signature = "metabCombiner",
    function(object, data = NULL)
{
      fdata <- featData(object, data = data)
      cols <- grep("^adduct", names(fdata), value = TRUE)
      return(fdata[c("rowID", cols)])
})

#' @rdname idData
#'
#' @export
setMethod("idData", signature = "metabCombiner", function(object, data = NULL)
{
    fdata <- featData(object, data = data)
    cols <- grep("^id", names(fdata), value = TRUE)
    return(fdata[c("rowID", cols)])
})

#' @rdname featData
#'
#' @export
setMethod("featData", signature = "metabCombiner", function(object, data = NULL)
{
    fdata <- object@featData
    if(is.null(data))
        return(fdata)
    else if(data == "x") data <- x(object)
    else if(data == "y") data <- y(object)
    else if(!(data %in% datasets(object)))
        stop("no dataset found with identifer in data argument")
    cols <- grep(paste("_",data,"$",sep = ""), names(fdata), value = TRUE)
    return(fdata[c("rowID", cols)])
})


#' @rdname combineData
#'
#' @export
setMethod("combineData", signature = "metabCombiner", function(object)
{
    fdata <- featData(object)
    samples_extras <- unlist(lapply(datasets(object), function(d)
        c(getSamples(object, d), getExtra(object, d))))
    seTable <- combinedTable(object)[samples_extras]
    table <- cbind.data.frame(fdata, seTable)
    return(table)
})


#' @rdname mzData
#'
#' @export
setMethod("mzData", signature = "metabCombiner",
    function(object, data = NULL, value = c("observed", "mean"))
{
    value <- match.arg(value)
    if(value == "mean") data <- NULL
    fdata <- featData(object, data = data)
    cols <- grep("^mz", names(fdata), value = TRUE)
    if(value == "observed")
        return(fdata[c("rowID", cols)])
    else
        return(matrixStats::rowMeans2(as.matrix(fdata[cols])))
})

#' @rdname QData
#'
#'@export
setMethod("QData", signature = "metabCombiner",
    function(object, data = NULL, value = c("observed", "mean"))
{
      value <- match.arg(value)
      if(value == "mean") data <- NULL
      fdata <- featData(object, data = data)
      cols <- grep("^Q", names(fdata), value = TRUE)
      if(value == "observed")
          return(fdata[c("rowID", cols)])
      else
          return(matrixStats::rowMeans2(as.matrix(fdata[cols])))
})

#' @rdname rtData
#'
#'@export
setMethod("rtData", signature = "metabCombiner",
    function(object, data = NULL, value = c("observed", "mean"))
{
    value <- match.arg(value)
    if(value == "mean") data <- NULL
    fdata <- featData(object, data = data)
    cols <- grep("^rt", names(fdata), value = TRUE)
    if(value == "observed")
        return(fdata[c("rowID", cols)])
    else
        return(matrixStats::rowMeans2(as.matrix(fdata[cols])))
})



