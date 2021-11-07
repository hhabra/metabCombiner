
#' @title Form a metabCombiner object.
#'
#' @description
#' This constructs an object of type \code{metabCombiner} from a pair of
#' metabolomics datasets, formatted as either \code{metabData} (single-dataset
#' class) or \code{metabCombiner} (combined-dataset class). An initial table of
#' possible feature pair alignments is constructed by grouping features into m/z
#' groups controlled by the binGap argument
#'
#' @param xdata metabData or metabCombiner object
#'
#' @param ydata metabData or metabCombiner object
#'
#' @param binGap numeric parameter used for grouping features by m/z.
#' See ?mzGroup for more details.
#'
#' @param xid character. If xdata is a \code{metabData}, assigns a new identifier
#' for this dataset; if xdata is a \code{metabCombiner}, selects one of the
#' existing dataset IDs to represent xdata. See details for more information.
#'
#' @param yid character. If ydata is a \code{metabData}, assigns a new identifier
#' for this dataset; if ydata is a \code{metabCombiner}, selects one of the
#' existing dataset IDs to represent ydata. See details for more information.
#'
#' @param means logical. Option to take average m/z, rt, and/or Q from
#' \code{metabComber}. May be a vector (length = 3), a single value (TRUE/FALSE),
#' or a list with names "mz", "rt", "Q" as names.
#'
#' @param rtOrder logical. If set to TRUE, retention order consistency expected
#' when resolving conflicting alignments for \code{metabCombiner} object inputs.
#'
#' @param impute logical. If TRUE, imputes the mean mz/rt/Q values for missing
#' features in \code{metabCombiner} object inputs before use in alignment (not
#' recommended for disparate data alignment); if FALSE, features with missing
#' information are dropped.
#'
#' @details
#' This function serves as a constructor of the \code{metabCombiner} combined
#' dataset class and the entry point to the main workflow for pairwise dataset
#' alignment. Two arguments must be specified, \code{xdata} and \code{ydata},
#' which must be either \code{metabData} or \code{metabCombiner} objects.
#' There are four scenarios listed here:
#'
#' 1) If xdata & ydata are \code{metabData} objects, a new \code{metabCombiner}
#' object is constructed with an alignment of this pair. New character
#' identifiers are assigned to each dataset (xid & yid, respectively); if these
#' are unassigned, then "1" and "2" will be their respective ids. xdata & ydata
#' will be the active "dataset x" and "dataset y" used for the paired alignment.
#'
#' 2) If xdata is a \code{metabCombiner} and ydata is a \code{metabData}, then
#' the result is the existing \code{metabCombiner} xdata augmented by an
#' additional dataset, ydata. One set of meta-data (id, m/z, rt, Q, adduct
#' labels) from xdata is used for alignment with the respective information
#' from ydata, which is controlled by the \code{xid} argument; see the
#' \code{\link{datasets}} method for extracting existing dataset ids. A new
#' identifier yid is assigned to ydata, which must be distinct from the current
#' dataset identifier.
#'
#' 3) If xdata is a \code{metabData} and ydata is a \code{metabCombiner}, then
#' a similar process to #2 occurs, with xdata augmented to the existing ydata
#' object and one of the constitutent dataset's meta-data is accessed, as
#' controlled by the yid argument. One major difference is that rts of ydata
#' serve as the "reference" or dependent variable in the spline-fitting step.
#'
#' 4) If xdata and ydata are both \code{metabCombiner} objects, the resulting
#' \code{metabCombiner} object aligns information from both combined datasets.
#' As before, one set of values contained in xdata (specified by xid argument)
#' is used to align to the values from ydata (controlled by yid argument).
#' The samples and extra columns are concatenated from all datasets.
#'
#' For \code{metabCombiner} object inputs, the full workflow
#' (\code{\link{selectAnchors}}, \code{\link{fit_gam}}/\code{\link{fit_loess}},
#' \code{\link{calcScores}}, \code{\link{labelRows}}) must be performed before
#' further alignment. If not completed already, features are pared down to
#' 1-1 alignments via the resolveConflicts approach (see: help(resolveRows)).
#' Features may not be used more than twice and will be removed if they are
#' detected as duplicates.
#'
#' The mean of the numeric fields (m/z, rt, Q) from all constituent datasets can
#' be used in alignment in place of values from a single dataset. These are
#' controlled by the means argument. By default this is a list value with "mz",
#' "rt" and "Q" as names, but may also accept a single logical or a length-3
#' logical vector. If set to a single logical value, then all three fields are
#' averaged (TRUE) or not averaged (FALSE). If a three-length argument is
#' supplied (e.g. c(TRUE, FALSE, FALSE)), then the values correspond to m/z, rt,
#' and Q respectively. RT averaging is generally not recommended for disparate
#' data alignment.
#'
#' If missing features have been incorporated into the \code{metabCombiner},
#' they an be imputed using the average m/z, rt, and Q values for that feature
#' in datasets in which it is present by setting \code{impute} to TRUE.
#' Likewise, this option is not recommended for disparate data alignment.
#'
#' @return a \code{metabCombiner} object constructed from xdata and ydata, with
#' features grouped by m/z according to the binGap argument.
#'
#' @note If using a \code{metabCombiner} object as input, only one row is
#' allowed per feature corresponding to its first appearance. It is strongly
#' recommended to reduce the table to 1-1 paired matches prior to aligning it
#' with a new dataset.
#'
#' @examples
#' data(plasma30)
#' data(plasma20)
#'
#' p30 <- metabData(plasma30, samples = "CHEAR")
#' p20 <- metabData(plasma20, samples = "Red", rtmax = 17.25)
#'
#' p.comb = metabCombiner(xdata = p30, ydata = p20, binGap = 0.0075,
#'                        xid = "p30", yid = "p20")
#'
#' @export
##
metabCombiner <- function(xdata, ydata, binGap = 0.005, xid = NULL, yid = NULL,
                        means = list('mz' = FALSE, 'rt' = FALSE, 'Q' = FALSE),
                        rtOrder = TRUE, impute = FALSE)
{
    check_combine_pars(binGap, means, xid, yid)
    if(methods::is(xdata,"metabData") & methods::is(ydata,"metabData"))
        object <- combine_default(xdata, ydata, binGap, xid, yid)
    else if(methods::is(xdata,"metabCombiner") & methods::is(ydata,"metabData"))
        object <- combine_X_y(xdata, ydata, binGap, xid, yid, means,
                              rtOrder, impute)
    else if(methods::is(xdata,"metabData") & methods::is(ydata,"metabCombiner"))
        object <- combine_x_Y(xdata, ydata, binGap, xid, yid, means,
                              rtOrder, impute)
    else if(methods::is(xdata,"metabCombiner") &
            methods::is(ydata,"metabCombiner"))
        object <- combine_X_Y(xdata, ydata, binGap, xid, yid, means,
                              rtOrder, impute)
    else
        stop("improper inputs for arguments 'xdata' & 'ydata'")

    return(object)
}

##xdata: metabData, ydata: metabData
combine_default <- function(xdata, ydata, binGap, xid, yid)
{
    combinerCheck(isMetabData(xdata), "metabData")
    combinerCheck(isMetabData(ydata), "metabData")
    object <- new("metabCombiner")
    xset <- getData(xdata)
    yset <- getData(ydata)
    xid <- ifelse(is.null(xid), "1", as.character(xid[1]))
    yid <- ifelse(is.null(yid), "2", as.character(yid[1]))
    if(xid %in% c("x", "y") | yid %in% c("x", "y"))
        stop("naming single datasets as 'x' or 'y' forbidden")
    xygroups <- mzGroup(xset = xset, yset = yset, binGap = binGap)
    samples <- list(); extra <- list();  nonmatched <- list()
    samples[[xid]] <- getSamples(xdata); samples[[yid]] <- getSamples(ydata)
    extra[[xid]] <- getExtra(xdata);     extra[[yid]] <- getExtra(ydata)
    nonmatched[[xid]] <- dplyr::filter(xygroups[["x"]],.data$group == 0)
    nonmatched[[yid]] <- dplyr::filter(xygroups[["y"]],.data$group == 0)
    xy <- list("x" = xid, "y" = yid, "xtype" = "metabData",
              "ytype" = "metabData")
    datasets <- list("x" = xid, "y" = yid)
    object <- update_mc(object, samples = samples, extra = extra, xy = xy,
                        stats = c("binGap","input_size_X", "input_size_Y"),
                        values = c(binGap, nrow(xset),nrow(yset)))
    xset <- dplyr::filter(xygroups[["x"]], .data$group > 0) %>%
            dplyr::arrange(.data$group)
    yset <- dplyr::filter(xygroups[["y"]], .data$group > 0) %>%
            dplyr::arrange(.data$group)
    nGroups <- max(xset[["group"]])
    object <- update_mc(object, nonmatched = nonmatched, datasets = datasets,
                    stats = c("grouped_X", "grouped_Y", "nGroups"),
                    values = c(nrow(xset), nrow(yset), nGroups))
    object <- form_tables(object, xset, yset, NULL, NULL, nGroups)
    return(object)
}

## xdata: metabCombiner,  ydata: metabData
combine_X_y <- function(xdata, ydata, binGap, xid, yid, means, rtOrder, impute)
{
    combinerCheck(isMetabCombiner(xdata), "metabCombiner")
    combinerCheck(isMetabData(ydata), "metabData")
    xid <- ifelse(is.null(xid), x(xdata), xid)
    if(xid == "x") xid <- x(xdata)
    else if(xid == "y") xid <- y(xdata)
    if(!(xid %in% datasets(xdata)))
        stop("no dataset with label ", xid, " found in xdata")
    xset <- form_dataset(xdata, data = xid, means = means, rtOrder = rtOrder,
                        impute = impute)
    yid <- ifelse(is.null(yid), as.character(length(datasets(xdata))+1), yid)
    if(yid %in% datasets(xdata))
        stop("argument 'yid' must be a character identifier not in xdata")
    if(yid %in% c("x", "y") )
        stop("naming single datasets as 'x' or 'y' forbidden")
    yset <- getData(ydata)
    xygroups <- mzGroup(xset = xset, yset = yset, binGap = binGap)
    samples <- getSamples(xdata);    samples[[yid]] <- getSamples(ydata)
    extra <- getExtra(xdata);        extra[[yid]] <- getExtra(ydata)
    nonmatched <- nonmatched(xdata, data = NULL)
    nonmatched[[yid]] <- dplyr::filter(xygroups[["y"]],.data$group == 0)
    xy <- list("x" = xid, "y" = yid, "xtype" = "metabCombiner",
               "ytype" = "metabData")
    datasets <- list("x" = datasets(xdata), "y" = yid)
    object <- methods::new("metabCombiner")
    object <- update_mc(object, samples = samples, extra = extra, xy = xy,
                        stats = c("binGap","input_size_X","input_size_Y"),
                        values = c(binGap, nrow(xset), nrow(yset)))
    xset <- dplyr::filter(xygroups[["x"]], .data$group > 0) %>%
            dplyr::arrange(.data$group)
    yset <- dplyr::filter(xygroups[["y"]], .data$group > 0) %>%
            dplyr::arrange(.data$group)
    nGroups <- max(xset[["group"]])
    object <- update_mc(object, nonmatched = nonmatched, datasets = datasets,
                        stats = c("grouped_X", "grouped_Y", "nGroups"),
                        values = c(nrow(xset), nrow(yset), nGroups))
    xfeat <- featdata(xdata)[match(xset[["rowID"]],featdata(xdata)[["rowID"]]),]
    object <- form_tables(object, xset, yset, xfeat, NULL, nGroups)
    return(object)
}

## xdata: metabData,  ydata: metabCombiner
combine_x_Y <- function(xdata, ydata, binGap, xid, yid, means, rtOrder, impute)
{
    combinerCheck(isMetabData(xdata), "metabData")
    combinerCheck(isMetabCombiner(ydata), "metabCombiner")
    xset <- getData(xdata)
    xid <- ifelse(is.null(xid), as.character(length(datasets(ydata))+1), xid)
    if(xid %in% datasets(ydata))
        stop("argument 'xid' must be a character identifier not in ydata")
    if(xid %in% c("x", "y"))
        stop("naming initial datasets as 'x' or 'y' forbidden")
    yid <- ifelse(is.null(yid), y(ydata), yid)
    if(yid == "x") yid <- x(ydata)
    if(yid == "y") yid <- y(ydata)
    if(!(yid %in% datasets(ydata)))
        stop("no dataset with label ", yid, " found in ydata")
    yset <- form_dataset(ydata, data = yid, means = means, rtOrder = rtOrder,
                        impute = impute)
    xygroups <- mzGroup(xset = xset, yset = yset, binGap = binGap)
    datasets = list("x" = xid, "y" = datasets(ydata))
    dats <- as.character(unlist(datasets))
    samples <- append(list(getSamples(xdata)), getSamples(ydata));
    extra <- append(list(getExtra(xdata)), getExtra(ydata));
    nn <- list(xid = c(dplyr::filter(xygroups[["x"]],.data$group == 0)))
    nn <- append(nn, nonmatched(ydata, data = NULL))
    names(samples) <- names(extra) <- names(nn) <- dats
    xy <- list("x" = xid, "y" = yid, "xtype" = "metabData",
               "ytype" = "metabCombiner")
    object <- methods::new("metabCombiner")
    object <- update_mc(object, samples = samples, extra = extra, xy = xy,
                        stats = c("binGap","input_size_X", "input_size_Y"),
                        values = c(binGap, nrow(xset),nrow(yset)))
    xset <- dplyr::filter(xygroups[["x"]], .data$group > 0) %>%
            dplyr::arrange(.data$group)
    yset <- dplyr::filter(xygroups[["y"]], .data$group > 0) %>%
            dplyr::arrange(.data$group)
    nGroups <- max(xset[["group"]])
    object <- update_mc(object, nonmatched = nn, datasets = datasets,
                        stats = c("grouped_X", "grouped_Y", "nGroups"),
                        values = c(nrow(xset), nrow(yset), nGroups))
    yfeat <- featdata(ydata)[match(yset[["rowID"]],featdata(ydata)[["rowID"]]),]
    object <- form_tables(object, xset, yset, NULL, yfeat, nGroups)
    return(object)
}

## xdata: metabCombiner, ydata: metabCombiner
combine_X_Y <- function(xdata, ydata, binGap, xid, yid, means, rtOrder, impute)
{
    combinerCheck(isMetabCombiner(xdata), "metabCombiner")
    combinerCheck(isMetabCombiner(ydata), "metabCombiner")
    xid <- ifelse(is.null(xid), x(xdata), xid)
    if(xid == "x") xid <- x(xdata)
    if(xid == "y") xid <- y(xdata)
    if(!(xid %in% datasets(xdata)))
        stop("no dataset with label ", xid, " found in xdata")
    xset <- form_dataset(xdata, data = xid, means = means, rtOrder = rtOrder,
                        impute = impute)
    yid <- ifelse(is.null(yid), y(ydata), yid)
    if(yid == "x") yid <- x(ydata)
    if(yid == "y") yid <- y(ydata)
    if(!(yid %in% datasets(ydata)))
        stop("no dataset with label ", yid, " found in ydata")
    yset <- form_dataset(ydata, data = yid, means = means, rtOrder = rtOrder,
                         impute = impute)
    xygroups <- mzGroup(xset = xset, yset = yset, binGap = binGap)
    dats <- c(datasets(xdata), datasets(ydata))
    if(any(duplicated(dats))){
        warning("duplicate dataset IDs detected; some original IDs modified")
        indices <- duplicated(dats)
        dats[indices] <- paste(dats[indices], 2, sep = ".")
    }
    samples <- c(getSamples(xdata), getSamples(ydata)); names(samples) <- dats
    extra <- c(getExtra(xdata), getExtra(ydata));  names(extra) <- dats
    nonmatched <- c(nonmatched(xdata, NULL), nonmatched(ydata, NULL))
    datasets <- list("x" = datasets(xdata),
                     "y" = dats[seq(length(datasets(xdata))+1,length(dats))])
    xy <- list("x" = xid, "y" = yid, "xtype" = "metabCombiner",
               "ytype" = "metabCombiner")
    object <- methods::new("metabCombiner")
    object <- update_mc(object, samples = samples, extra = extra, xy = xy,
                        stats = c("binGap","input_size_X", "input_size_Y"),
                        values = c(binGap, nrow(xset),nrow(yset)))
    xset <- dplyr::filter(xygroups[["x"]], .data$group > 0) %>%
            dplyr::arrange(.data$group)
    yset <- dplyr::filter(xygroups[["y"]], .data$group > 0) %>%
            dplyr::arrange(.data$group)
    nGroups <- max(xset[["group"]])
    object <- update_mc(object, nonmatched = nonmatched, datasets = datasets,
                        stats = c("grouped_X", "grouped_Y", "nGroups"),
                        values = c(nrow(xset), nrow(yset), nGroups))
    xfeat <- featdata(xdata)[match(xset[["rowID"]],featdata(xdata)[["rowID"]]),]
    yfeat <- featdata(ydata)[match(yset[["rowID"]],featdata(ydata)[["rowID"]]),]
    object <- form_tables(object, xset, yset, xfeat, yfeat, nGroups)
    return(object)
}
