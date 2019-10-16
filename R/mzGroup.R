##m/z grouping function

#' @title Binning of mass spectral features in m/z dimension
#'
#' @description
#' Features in two input feature lists are grouped by their m/z values.
#'
#' @details
#' The m/z values from both datasets are pooled, sorted, and binned by the
#' \code{binGap} argument. Feature groups form when there is at least one pair
#' of features from both datasets whose consecutive difference is less than
#' \code{binGap}.
#'
#' @param object metabCombiner object
#'
#' @param binGap numeric. Gap between concecutive sorted pooled feature m/z values.
#'
#' @return metabCombiner object with group information in data objects
##
mzGroup <- function(object, binGap = 0.005){
    if(class(object) != "metabCombiner")
        stop(base::paste(object, "is not a metabCombiner object"), sep = " ")

    xdata = dplyr::mutate(getData(object, "x"), index = 1:nrow(getData(object, "x")))
    ydata = dplyr::mutate(getData(object, "y"), index = 1:nrow(getData(object, "y")))

    if(nrow(xdata) == 0 | nrow(ydata) == 0)
        stop(base::paste("missing xdata or ydata in object", object, sep = " "))

    if(!("mz" %in% names(xdata)) | !("mz" %in% names(ydata)))
        stop("m/z column missing in either xdata or ydata")

    if(is.null(binGap) | length(binGap) == 0)
        stop("parameter 'binGap' must be defined")

    if(binGap <= 0 | !is.numeric(binGap))
        stop("parameter 'binGap' must be a positive numeric constant")

    if(binGap >= 0.3)
        stop("parameter 'binGap' is too high")

    if(binGap > 0.05)
        warning("parameter 'binGap' is very high. binGap < 0.05 is recommended.")

    object@binGap = binGap

    #forming m/z groups table
    mzGroups <- rbind(dplyr::select(xdata, mz, index),
                      dplyr::select(ydata, mz, index)) %>%
                dplyr::mutate(set = c(rep("x", nrow(xdata)),
                                      rep("y", nrow(ydata)))) %>%
                dplyr::arrange(mz)

    if(!is.numeric(mzGroups[["mz"]]) | any(is.na(mzGroups[["mz"]])) |
        any(mzGroups[["mz"]] < 0))
        stop("At least one negative, missing, or non-numeric m/z value")

    mzGroups[["groups"]] <- .Call("binByMZ",
                                  mz = mzGroups[["mz"]],
                                  datasets = mzGroups[["set"]],
                                  gap = binGap,
                                  PACKAGE = "Combiner")

    mzGroupsX <- dplyr::filter(mzGroups, set == "x") %>% dplyr::arrange(index)
    mzGroupsY <- dplyr::filter(mzGroups, set == "y") %>% dplyr::arrange(index)

    object@xdata@data[["group"]] <- mzGroupsX[["groups"]]
    object@ydata@data[["group"]] <- mzGroupsY[["groups"]]

    return(object)
}


