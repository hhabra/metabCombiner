##m/z grouping function

#' @title Binning of mass spectral features in m/z dimension
#'
#' @description
#' Features in two input feature lists are grouped by their m/z values.
#'
#' @param xset data frame containing metabolomics features
#'
#' @param yset data frame containing metabolomics features
#'
#' @param binGap numeric gap value between consecutive sorted & pooled feature
#'              m/z values.
#'
#' @details
#' The m/z values from both datasets are pooled, sorted, and binned by the
#' \code{binGap} argument. Feature groups form when there is at least one pair
#' of features from both datasets whose consecutive difference is less than
#' \code{binGap}. Grouped features are joined together in \code{combinedTable}
#' data report.
#'
#' @return list object containing updated xset & yset with group information
mzGroup <- function(xset, yset, binGap){
    #forming m/z groups table
    mzGroups <- data.frame(`mz` = c(xset[["mz"]], yset[["mz"]])) %>%
                dplyr::mutate(`set` = c(rep("x", nrow(xset)),
                                        rep("y", nrow(yset))),
                            `index` = c(seq(1,nrow(xset)),
                                        seq(1,nrow(yset)))) %>%
                dplyr::arrange(.data$mz)

    if(!is.numeric(mzGroups[["mz"]]) | any(is.na(mzGroups[["mz"]])) |
        any(mzGroups[["mz"]] < 0))
        stop("At least one negative, missing, or non-numeric m/z value")

    mzGroups[["group"]] <- .Call("binByMZ",
                                    mz = mzGroups[["mz"]],
                                    datasets = mzGroups[["set"]],
                                    gap = binGap,
                                    PACKAGE = "metabCombiner")

    xgroups <- dplyr::filter(mzGroups, .data$set == "x") %>%
                dplyr::arrange(.data$index)
    ygroups <- dplyr::filter(mzGroups, .data$set == "y") %>%
                dplyr::arrange(.data$index)

    xset[["group"]] <- xgroups[["group"]]
    yset[["group"]] <- ygroups[["group"]]

    return(list(x = xset, y = yset))
}


