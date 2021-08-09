##Methods for extracting feature meta-data
#' @include generics.R classes.R


#' @rdname adductdata
#'
#' @export
setMethod("adductdata",signature = "metabCombiner",
    function(object, data = NULL)
{
    fdata <- featdata(object)
    if(is.null(data))
        data <- datasets(object)
    else if(data == "x") data <- x(object)
    else if(data == "y") data <- y(object)
    else if(!(data %in% datasets(object)))
        stop("no dataset identifer found for 'data' argument")
    cols <- paste("adduct", data, sep = "_")
    return(fdata[c("rowID", cols)])
})

#' @rdname iddata
#'
#' @export
setMethod("iddata", signature = "metabCombiner", function(object, data = NULL)
{
    fdata <- featdata(object)
    if(is.null(data))
        data <- datasets(object)
    else if(data == "x") data <- x(object)
    else if(data == "y") data <- y(object)
    else if(!(data %in% datasets(object)))
        stop("no dataset found with identifer in data argument")
    cols <- paste("id", data, sep = "_")
    return(fdata[c("rowID", cols)])
})

#' @rdname featdata
#'
#' @export
setMethod("featdata", signature = "metabCombiner", function(object, data = NULL)
{
    fdata <- object@featdata
    if(is.null(data))
        return(fdata)
    else if(data == "x") data <- x(object)
    else if(data == "y") data <- y(object)
    else if(!(data %in% datasets(object)))
        stop("no dataset found with identifer in data argument")
    cols <- grep(paste("_",data,"$",sep = ""), names(fdata), value = TRUE)
    return(fdata[c("rowID", cols)])
})


#' @rdname mzdata
#'
#' @export
setMethod("mzdata", signature = "metabCombiner",
    function(object, data = NULL, value = c("obs", "mean"))
{
    value = match.arg(value)
    fdata <- featdata(object)
    if(is.null(data) | value == "mean") data <- datasets(object)
    else if(data == "x") data <- x(object)
    else if(data == "y") data <- y(object)
    else if(!(data %in% datasets(object)))
        stop("no dataset identifer found for 'data' argument")

    cols <- paste("mz", data, sep = "_")
    if(value == "obs")
        return(fdata[c("rowID", cols)])
    else
        return(matrixStats::rowMeans2(as.matrix(fdata[cols])))
})

#' @rdname Qdata
#'
#'@export
setMethod("Qdata", signature = "metabCombiner",
    function(object, data = NULL, value = c("obs", "mean"))
{
    value = match.arg(value)
    fdata <- featdata(object)
    if(is.null(data) | value == "mean")
        data <- datasets(object)
    else if(data == "x") data <- x(object)
    else if(data == "y") data <- y(object)
    else if(!(data %in% datasets(object)))
        stop("no dataset identifer found for 'data' argument")
    cols <- paste("Q", data, sep = "_")

    if(value == "obs")
        return(fdata[c("rowID", cols)])
    else
        return(matrixStats::rowMeans2(as.matrix(fdata[cols])))
})

#' @rdname rtdata
#'
#'@export
setMethod("rtdata", signature = "metabCombiner",
    function(object, data = NULL, value = c("obs", "mean"))
{
    value = match.arg(value)
    fdata <- featdata(object)
    if(is.null(data) | value == "mean")
        data <- datasets(object)
    else if(data == "x") data <- x(object)
    else if(data == "y") data <- y(object)
    else if(!(data %in% datasets(object)))
        stop("no dataset identifer found for 'data' argument")
    cols <- paste("rt", data, sep = "_")

    if(value == "obs")
        return(fdata[c("rowID", cols)])
    else
        return(matrixStats::rowMeans2(as.matrix(fdata[cols])))
})



