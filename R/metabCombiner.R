##
#' @title Form a metabCombiner object.
#'
#' @description
#' \code{metabCombiner()} takes two metabolomics featurelists, contained in
#' \code{metabData} objects,and constructs a merged dataset containing groups of
#' features detected in both featurelists with similar m/z values.
#'
#' @param xdata metabData object. One of two datasets to be Combined.
#'
#' @param ydata metabData object. One of two datasets to be Combined.
#'
#' @param binGap numeric. Parameter used for grouping features by m/z.
#' See ?mzGroup for more details.
#'
#' @details
#' The \code{binGap} argument defines consecutive differences between grouped
#' m/z features. \code{metabCombiner} objects are used for all subsequent
#' processing steps.
#'
#' @return a metabCombiner object
#'
#' @examples
#' \dontrun{
#' library(metabCombiner)
#' data(plasma30)
#' data(plasma20)
#'
#' p30 <- metabData(plasma30, samples = "CHEAR")
#' p20 <- metabData(plasma20, samples = "Red", rtmax = 17.25)
#'
#' p.combined = metabCombiner(xdata = p30, ydata = p20, binGap = 0.0075)
#' }
#' @export
##
metabCombiner <- function(xdata, ydata, binGap = 0.005){
    codeX = isMetabData(xdata)
    codeY = isMetabData(ydata)

    if(codeX)
        stop(paste("invalid xdata parameter:",combinerError(codeX, "metabData")))

    if(codeY)
        stop(paste("invalid ydata parameter:",combinerError(codeY, "metabData")))

    if(is.null(binGap) | length(binGap) == 0)
        stop("parameter 'binGap' must be defined")

    if(binGap <= 0 | !is.numeric(binGap))
        stop("parameter 'binGap' must be a positive numeric constant")

    if(binGap >= 0.25)
        stop("parameter 'binGap' is too large")

    if(binGap > 0.05)
        warning("parameter 'binGap' is very high. binGap < 0.05 is recommended.")

    object <- new("metabCombiner")
    object@stats[["binGap"]] = binGap

    xset = getData(xdata)
    yset = getData(ydata)

    xygroups = mzGroup(xset = xset,
                       yset = yset,
                       binGap = binGap)

    object@stats[["input_size_X"]] = nrow(xset)
    object@nonmatched[["x"]] = dplyr::filter(xygroups[["x"]], .data$group == 0)
    xset = dplyr::filter(xygroups[["x"]], .data$group > 0) %>%
           dplyr::arrange(.data$group)
    object@stats[["grouped_size_X"]] = nrow(xset)

    object@stats[["input_size_Y"]] = nrow(yset)
    object@nonmatched[["y"]] = dplyr::filter(xygroups[["y"]], .data$group == 0)
    yset = dplyr::filter(xygroups[["y"]], .data$group > 0) %>%
           dplyr::arrange(.data$group)
    object@stats[["grouped_size_Y"]] = nrow(yset)

    nGroups = max(xset[["group"]])
    object@stats[["nGroups"]] = nGroups

    object@samples[["x"]] = xdata@samples
    object@samples[["y"]] = ydata@samples

    object@extra[["x"]] = xdata@extra
    object@extra[["y"]] = ydata@extra

    object = formCombinedTable(object, xset, yset, nGroups)

    return(object)
}

