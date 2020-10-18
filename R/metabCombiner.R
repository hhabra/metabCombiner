##
#' @title Form a metabCombiner object.
#'
#' @description
#' \code{metabCombiner()} takes two metabolomics featurelists, contained in
#' \code{metabData} objects and constructs a merged dataset containing groups
#' of features detected in both featurelists with similar m/z values.
#'
#' @param xdata metabData object. One of two datasets to be combined.
#'
#' @param ydata metabData object. One of two datasets to be combined.
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
#' data(plasma30)
#' data(plasma20)
#'
#' p30 <- metabData(plasma30, samples = "CHEAR")
#' p20 <- metabData(plasma20, samples = "Red", rtmax = 17.25)
#'
#' p.comb = metabCombiner(xdata = p30, ydata = p20, binGap = 0.0075)
#'
#' @export
##
metabCombiner <- function(xdata, ydata, binGap = 0.005){
    combinerCheck(isMetabData(xdata), "metabData")
    combinerCheck(isMetabData(ydata), "metabData")

    if(is.null(binGap) | length(binGap) == 0)
        stop("argument 'binGap' must be defined")

    if(binGap <= 0 | !is.numeric(binGap))
        stop("argument 'binGap' must be a positive numeric constant")

    if(binGap >= 0.1)
        stop("argument 'binGap' is too large")

    if(binGap > 0.05)
        warning("large 'binGap' argument value; binGap < 0.01 is recommended")

    object <- new("metabCombiner")
    xset <- getData(xdata)
    yset <- getData(ydata)
    xygroups <- mzGroup(xset = xset, yset = yset, binGap = binGap)

    object <- update_mc(object, data = "x", samples = getSamples(xdata),
                    extra = getExtra(xdata),
                    stats = c("binGap","input_size_X", "input_size_Y"),
                    values = c(binGap, nrow(xset),nrow(yset)),
                    nonmatched = dplyr::filter(xygroups$x,.data$group == 0))
    xset <- dplyr::filter(xygroups[["x"]], .data$group > 0) %>%
            dplyr::arrange(.data$group)
    yset <- dplyr::filter(xygroups[["y"]], .data$group > 0) %>%
            dplyr::arrange(.data$group)
    nGroups <- max(xset[["group"]])
    object <- update_mc(object, data = "y", samples = getSamples(ydata),
                    extra = getExtra(ydata),
                    stats = c("grouped_size_X", "grouped_size_Y", "nGroups"),
                    values = c(nrow(xset), nrow(yset), nGroups))

    object <- formCombinedTable(object, xset, yset, nGroups)
    return(object)
}

