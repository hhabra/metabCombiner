##Functions for selecting ordered pairs for retention time fitting

##
#' @title Select Matching Ids as Anchors
#'
#' @description
#' This is an optional helper function for \code{selectAnchors}. Uses identities
#' to guide selection of ordered retention time pairs. If useID option is set to
#' TRUE, it will select pairs of features with matching ID character strings
#' before proceeding with iterative anchor selection.
#'
#' @param cTable  data frame, contains only feature ids, mzs, rts, Qs, & labels
#'
#' @param windx   numeric positive retention time exclusion window in X dataset.
#'
#' @param windy   numeric positive retention time exclusion windown in Y dataset.
#'
#' @param useID   logical. Operation proceeds if TRUE, terminates otherwise.
#'
#' @param brackets If useID = TRUE, bracketed identity strings of the types
#' included in this argument will be ignored
#'
#' @note Identity anchors are allowed to violate constraints of m/z, Q, and rtq
#' difference tolerances, and will not be removed if they fall within the rt
#' exclusion window of other features. If a name appears more than once, only
#' the pair with the highest relative abundance is selected.
#'
#' @seealso
#' \code{\link{selectAnchors}}
#'
##
identityAnchorSelection <- function(cTable, windx, windy, useID, brackets)
{
    if(!useID)
        return(cTable)

    cTable[["labels"]] = compare_strings(cTable[["idx"]],
                                         cTable[["idy"]],
                                         "I", "P", brackets)

    if(any(cTable[["labels"]] == "I"))
    {
        #duplicate identities removed; most abundant representation kept
        xTable = dplyr::mutate(cTable,
                               `row` = 1:nrow(cTable),
                               `idx` = tolower(.data$idx)) %>%
                 dplyr::filter(.data$labels == "I") %>%
                 dplyr::arrange(dplyr::desc(.data$Qx+ .data$Qy)) %>%
                 dplyr::filter(duplicated(.data$idx))

        cTable$labels[xTable[["row"]]] = "N"

        ids = which(cTable[["labels"]] == "I")

        cTable[["labels"]] = .Call("selectAnchorsByID",
                                    labels = cTable[["labels"]],
                                    ids = ids,
                                    rtx = cTable[["rtx"]],
                                    rty = cTable[["rty"]],
                                    windx = windx,
                                    windy = windy,
                                    PACKAGE = "metabCombiner")
    }

    return(cTable)

}

#' @title Iteratively Select Ordered Feature Pairs
#'
#' @description
#' This is a helper function for \code{selectAnchors}. Anchors are iteratively
#' selected from highly abundant feature pairs, subject to feature m/z, rt, & Q
#' constraints set by the user.
#'
#' @param cTable  data frame, contains only feature ids, mzs, rts, Qs, & labels
#'
#' @param windx   numeric positive retention time exclusion window in X dataset.
#'
#' @param windy   numeric positive retention time exclusion windown in Y dataset.
#'
#' @param swap  logical. When FALSE, searches for abundant features in dataset X,
#' complemented by dataset Y features; when TRUE, searches for abundant features
#' in dataset Y, complemented by dataset X features.
#'
#' @return
#' data frame of anchor feature alignments.
#'
#'@seealso
#' \code{\link{selectAnchors}}
##
iterativeAnchorSelection <- function(cTable, windx, windy, swap = FALSE){
    if(swap)
        cTable[["labels"]] = .Call("selectIterativeAnchors",
                                    labels = cTable[["labels"]],
                                    rtx = cTable[["rty"]],
                                    rty = cTable[["rtx"]],
                                    Qx = cTable[["Qy"]],
                                    Qy = cTable[["Qx"]],
                                    windx = windy,
                                    windy = windx,
                                    PACKAGE = "metabCombiner")
    else
        cTable[["labels"]] = .Call("selectIterativeAnchors",
                                    labels = cTable[["labels"]],
                                    rtx = cTable[["rtx"]],
                                    rty = cTable[["rty"]],
                                    Qx = cTable[["Qx"]],
                                    Qy = cTable[["Qy"]],
                                    windx = windx,
                                    windy = windy,
                                    PACKAGE = "metabCombiner")


    cTable <- dplyr::filter(cTable, .data$labels != "N")

    return(cTable)
}

#' @title Select Anchors for Nonlinear RT Model
#'
#' @description
#' A subset of all possible alignments in the \code{combinedTable} are used as
#' ordered pairs to anchor a retention time projection model. Alignments of
#' abundant features are prominent targets for anchor selection, but shared
#' identified features (i.e. feature pairs where idx = idy) may also be used.
#'
#' @param object metabCombiner object.
#'
#' @param useID logical. Option to first search for IDs as anchors.
#'
#' @param tolmz numeric. m/z tolerance for prospective anchors
#'
#' @param tolQ numeric. Quantile Q tolerance for prospective anchors
#'
#' @param tolrtq numeric. Linear RT quantile tolerance for prosepctive anchors.
#'
#' @param windx numeric. Retention time exclusion window around each anchor in
#' X dataset. Optimal values are between 0.01 and 0.05 min (1-3s)
#'
#' @param windy numeric. Retention time exclusion window around each anchor in
#' dataset Y. Optimal values are between 0.01 and 0.05 min (1-3s)
#'
#' @param brackets_ignore If useID = TRUE, bracketed identity strings of the types
#' included in this argument will be ignored.
#'
#' @details
#' In order to map between two sets of retention times, a set of ordered pairs
#' need to be selected for the spline fit. This function relies on mutually
#' abundant features to select these ordered pairs. In iterative steps, the most
#' abundant (as indicated by Q value) in one dataset is selected along with its
#' counterpart, and all features within some retention time window specified by
#' \code{windx} & \code{windy} arguments are excluded. This process is repeated
#' until all features have been considered.
#'
#' \code{tolQ} & \code{tolmz} arguments restrict to feature pairs that have
#' differences in Q & m/z within these tolerances. \code{tolrtq} further limits
#' feature pairs to those with absolute differences in linear retention time
#' quantiles, here calculated as rtqx = (rtx - min(rtx)) / (max(rtx) - min(rtx))
#' and rtqy = (rty - min(rty)) / (max(rty) - min(rty))
#'
#' Shared identities (in which idx & idy columns have matching, non-empty &
#' non-bracketed strings) may be used if \code{useID} ise set to \code{TRUE}. In
#' this case, shared identities will be searched first and will not be subject
#' to any of the restrictions in m/z, Q, or rt. The iterative process proceeds
#' after processing of shared identities.
#'
#' @return \code{metabCombiner} object with updated \code{anchors} slot. This
#' is a data.frame of feature pairs that shall be used to map between retention
#' times using a nonlinear spline model.
#'
#' \item{idx}{identities of features from dataset X}
#' \item{idy}{identities of features from dataset Y}
#' \item{mzx}{m/z values of features from dataset X}
#' \item{mzy}{m/z values of features from dataset Y}
#' \item{rtx}{retention time values of features from dataset X}
#' \item{rty}{retention time values of features from dataset Y}
#' \item{Qx}{abundance quantile values of features from dataset X}
#' \item{Qy}{abundance quantile values of features from dataset Y}
#' \item{adductX}{adduct label of features from dataset X}
#' \item{adductY}{adduct label of features from dataset Y}
#' \item{group}{m/z feature group of feature pairing}
#' \item{labels}{anchor labels; "I" for identity, "A" for normal anchors}
#'
#' @examples
#' \dontrun{
#' library(metabCombiner)
#' data(plasma30)
#' data(plasma20)
#'
#' p30 <- metabData(plasma30, samples = "CHEAR")
#' p20 <- metabData(plasma20, samples = "Red", rtmax = 17.25)
#' p.combined = metabCombiner(xdata = p30, ydata = p20, binGap = 0.0075)
#'
#' ##example 1 (no known IDs used)
#' p.combined = selectAnchors(p.combined, tolmz = 0.003, tolQ = 0.3, windx = 0.03,
#'                            windy = 0.02)
#'
#' ##example 2 (known IDs used)
#' p.combined = selectAnchors(p.combined, useID = TRUE, tolmz = 0.003, tolQ = 0.3)
#'
#' ##To View Plot of Ordered Pairs
#' anchors = getAnchors(p.combined)
#' plot(anchors$rtx, anchors$rty, main = "Selected Anchor Ordered Pairs",
#'      xlab = "rtx", ylab = "rty")
#' }
#'
#' @seealso
#' \code{\link{getAnchors}}, \code{\link{fit_gam}}, \code{\link{fit_loess}}
#'
#' @export
selectAnchors <- function(object, useID = FALSE, tolmz = 0.003, tolQ = 0.3,
                          tolrtq = 0.5, windx = 0.03, windy = 0.03,
                          brackets_ignore = c("(", "[", "{"))
{
    code = isMetabCombiner(object)

    if(code)
        stop(combinerError(code, "metabCombiner"))

    if(!is.numeric(tolmz) | !is.numeric(tolQ)| !is.numeric(tolrtq) |
       !is.numeric(windx) | !is.numeric(windy))
        stop("expected numeric constants for tolmz, tolQ, tolrtq, windx, windy")

    if(tolmz <= 0 | tolQ <= 0 | tolrtq <= 0 | windx <= 0 | windy <= 0)
        stop("tolmz, tolQ, tolrtq, windx, windy must be positive constants")

    if(!is.logical(useID))
        stop("expected logical value for argument 'useID'")

    cTable <- combinedTable(object)[,1:15]

    #retention time extrema
    rtex = c(min(cTable[["rtx"]]), min(cTable[["rty"]]),
             max(cTable[["rtx"]]), max(cTable[["rty"]]))

    cTable = dplyr::select(cTable,
                           -.data$score,
                           -.data$rankX,
                           -.data$rankY,
                           -.data$rtProj) %>%
             dplyr::mutate(`labels` = rep("P", nrow(cTable)),
                           `rtqx` = (.data$rtx - rtex[1])/ (rtex[3] - rtex[1]),
                           `rtqy` = (.data$rty - rtex[2])/ (rtex[4] - rtex[2]))

    cTable = identityAnchorSelection(cTable,
                                     windx = windx,
                                     windy = windy,
                                     useID = useID,
                                     brackets = brackets_ignore)

    cTable = dplyr::filter(cTable, (abs(.data$mzx - .data$mzy) < tolmz &
                                    abs(.data$Qx - .data$Qy) < tolQ &
                                    abs(.data$rtqx - .data$rtqy) < tolrtq &
                                    .data$labels != "N") |
                                    .data$labels == "I")

    anchorlistXY <- iterativeAnchorSelection(cTable = cTable,
                                             windx = windx,
                                             windy = windy,
                                             swap = FALSE)
    anchorlistYX <- iterativeAnchorSelection(cTable = cTable,
                                             windx = windx,
                                             windy = windy,
                                             swap = TRUE)

    anchorlist <- dplyr::intersect(anchorlistXY, anchorlistYX)

    if(nrow(anchorlist) < 20)
        warning("Number of anchors is less than 20. Consider using
                different parameters or using identities.")

    object@anchors = anchorlist

    return(object)
}



