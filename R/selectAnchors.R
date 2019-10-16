##Functions for selecting ordered pairs for retention time fitting.

##
#' @title Combine unequal names or return identical strings
#'
#' @param s1  The first of two strings
#'
#' @param s2  The second of two strings
#'
#' @return semicolon-separated string (if s1 =/= s2); s1 (if s1 matches s2)
#'
#' @noRd
##
combine_strings <- function(s1, s2){
    s1 = ifelse(is.na(s1), "", s1)
    s2 = ifelse(is.na(s2), "", s2)

    ifelse(tolower(s1) == tolower(s2), s1, paste(s1,s2, sep = ";"))
}

##
#' @title Select Matching Ids as Anchors
#'
#' @description
#' This is an optional helper function for \code{selectAnchors}. Uses identities
#' to guide selection of ordered retention time pairs. If useID option is set to
#' true, it will select pairs of features with matching ID character strings before
#' proceeding with iterative selection.
#'
#' @param cTable  data frame, contains only feature ids, mzs, rts, Qs, & labels
#'
#' @param windX   numeric positive retention time exclusion window in X dataset.
#'
#' @param windY   numeric positive retention time exclusion windown in Y dataset.
#'
#' @param useID   logical. Operation proceeds if TRUE, terminates otherwise.
#'
#' @note Identity anchors may violate constraints of m/z and Q difference
#' tolerances, and will not be removed if they fall within the exclusion window
#' of other features. If a name appears more than once, only one pair will be
#' selected, and the others excluded.
#'
#' @seealso
#' [selectAnchors], [selectIterativeAnchors]
#'
##
identityAnchorSelection <- function(cTable, windX, windY, useID){
    if(!useID)
        return(cTable)

    length = nrow(cTable)

    ids = combine_strings(cTable[["idx"]], cTable[["idy"]])

    #mismatched and empty IDs
    mismatchIds = grep(";", ids)
    emptyIds = which(ids == "")
    nonIds = c(emptyIds, mismatchIds)

    matches = which(!((1:length) %in% nonIds))
    duplicateIds = matches[duplicated(tolower(cTable$idx[matches]))]

    if(length(matches) > 0){
        cTable$labels[matches] = "I"

        cTable$labels[duplicateIds] = "N"

        cTable[["labels"]] = .Call("selectAnchorsByID",
                                    labels = cTable[["labels"]],
                                    Ids = matches,
                                    rtx = cTable[["rtx"]],
                                    rty = cTable[["rty"]],
                                    windX = windX,
                                    windY = windY,
                                    PACKAGE = "Combiner")
    }

    return(cTable)

}

#' @title Iteratively Select Ordered Feature Pairs
#'
#' @description
#'
#'
#' @param cTable  data frame, contains only feature ids, mzs, rts, Qs, & labels
#'
#' @param windX   numeric positive retention time exclusion window in X dataset.
#'
#' @param windY   numeric positive retention time exclusion windown in Y dataset.
#'
#' @param useID   logical. Operation proceeds if TRUE, terminates early otherwise.
#'
#' @return
#' Data frame
#'
##
iterativeAnchorSelection <- function(cTable, windX, windY, swap = FALSE){
    if(swap)
        cTable[["labels"]] = .Call("selectIterativeAnchors",
                                    labels = cTable[["labels"]],
                                    rtx = cTable[["rty"]],
                                    rty = cTable[["rtx"]],
                                    Qx = cTable[["Qy"]],
                                    Qy = cTable[["Qx"]],
                                    windX = windY,
                                    windY = windX,
                                    PACKAGE = "Combiner")
    else
        cTable[["labels"]] = .Call("selectIterativeAnchors",
                                    labels = cTable[["labels"]],
                                    rtx = cTable[["rtx"]],
                                    rty = cTable[["rty"]],
                                    Qx = cTable[["Qx"]],
                                    Qy = cTable[["Qy"]],
                                    windX = windX,
                                    windY = windY,
                                    PACKAGE = "Combiner")


    cTable <- dplyr::filter(cTable, labels != "N")

    return(cTable)
}


#' @title Select Ordered Pair Features for Nonlinear RT Model
#'
#' @param object metabCombiner object.
#'
#' @param useID logical. Option to first search for IDs as anchors.
#'
#' @param tolMZ numeric. m/z tolerance for prospective anchors.
#'
#' @param tolQ numeric. Quantile Q tolerance for prospective anchors.
#'
#' @param windX numeric. Retention time exclusion window around each anchor in
#' X dataset.
#'
#' @param windY numeric. Retention time exclusion window around each anchor in
#' Y dataset.
#'
#' @return an updated metabCombiner object containing a table of ordered feature
#' pairs for fitting a nonlinear rt fitting model.
#'
#' @examples
#'
#'
#' @export
selectAnchors <- function(object, useID = FALSE, tolMZ = 0.005, tolQ = 0.5,
                          windX = 0.03, windY = 0.03){

    code = isMetabCombiner(object)

    if(code)
        stop(combinerError(code, "metabCombiner"))

    if(tolMZ <= 0 | tolQ <= 0 | windX <= 0 | windY <= 0)
        stop("Parameters tolMZ, tolQ, windX, windY must be positive constants")

    if(!is.logical(useID))
        stop("Expected logical value for argument 'useID'")

    cTable <- combinerTable(object)
    cTable = cTable[,1:15]

    cTable = dplyr::select(cTable, idx, idy, mzx, mzy,rtx, rty,
                           Qx, Qy, adductx, adducty,group) %>%
             dplyr::mutate(mzdiff = abs(mzx - mzy),
                           Qdiff = abs(Qx - Qy),
                           labels = rep("P", nrow(cTable))) %>%
             identityAnchorSelection(windX = windX,
                                     windY = windY,
                                     useID = useID) %>%
             dplyr::filter((abs(mzdiff) < tolMZ & abs(Qdiff) < tolQ) |
                             labels == "I") %>%
             dplyr::filter(labels != "N") %>%
             dplyr::select(-c(mzdiff, Qdiff))


    anchorlistXY <- iterativeAnchorSelection(cTable = cTable,
                                             windX = windX,
                                             windY = windY,
                                             swap = FALSE)
    anchorlistYX <- iterativeAnchorSelection(cTable = cTable,
                                             windX = windX,
                                             windY = windY,
                                             swap = TRUE)

    anchorlist <- dplyr::intersect(anchorlistXY, anchorlistYX)

    if(nrow(anchorlist) < 20)
        warning("Number of anchors is less than 20. Consider using
                different parameters or using identities.")

    object@anchors = anchorlist

    return(object)
}





