
#' Helper function containing C function call
#'
#' @noRd
get_labels <- function(fields, minScore, delta, maxRankX, maxRankY, method,
                    rterr, maxRTerr, balanced)
{
    fields[["labels"]] <- .Call("labelRows",
                                labels = fields$labels,
                                subgroup = fields$subgroup,
                                alt = fields$alt, mzx = fields$mzx,
                                mzy = fields$mzy, rtx = fields$rtx,
                                rty = fields$rty, score = fields$score,
                                rankX = fields$rankX, rankY = fields$rankY,
                                group = fields$group, balanced = balanced,
                                delta = delta, minScore = minScore,
                                maxRankX = as.integer(maxRankX[1]),
                                maxRankY = as.integer(maxRankY[1]),
                                method = method,maxRTerr = as.numeric(maxRTerr),
                                rterr = rterr, PACKAGE = "metabCombiner")

    return(fields)
}


#' Separate combinedTable into fields & values
#'
#' @noRd
prepare_fields_values <- function(cTable, useID, brackets_ignore)
{
    cTable <- cTable[with(cTable, order(`group`, desc(`score`))), ]
    values <- cTable[,-seq(1,length(combinerNames()))]
    fields <- cTable[combinerNames()]
    if(useID == TRUE)
        fields[["labels"]] <- compare_strings(fields[["idx"]], fields[["idy"]],
                                              "IDENTITY", "", brackets_ignore)
    else
        fields[["labels"]] <- ""
    fields[["subgroup"]] <- integer(nrow(fields))
    fields[["alt"]] <- integer(nrow(fields))

    if(any(cTable[["group"]]) == 0)
        group0 <- list(fields = fields[cTable[["group"]] == 0,],
                       values = values[cTable[["group"]] == 0,])
    else
        group0 <- NULL

    fields <- fields[cTable[["group"]] > 0,]
    values <- values[cTable[["group"]] > 0,]

    return(list(fields = fields, values = values, group0 = group0))
}


#' @title Annotate and Remove Report Rows
#'
#' @description
#' This is a method for annotating removable, conflicting, and identity-matched
#' feature pair alignment (FPA) rows in the \code{combinedTable} report. Simple
#' thresholds for score, rank, retention time error and delta score can
#' computationally reduce the set of possible FPAs to the most likely feature
#' matches. FPAs falling within some small delta score or mz/rt of the
#' top-ranked pair are organized into subgroups to facilitate inspection.
#' Automated reduction to 1-1 pairs is also possible with this function.
#'
#' \code{reduceTable} behaves identically to labelRows, but with a focus on
#' automated table reduction. Rank threshold defaults in \code{reduceTable} are
#' also stricter than in \code{labelRows}.
#'
#' @param object Either a \code{metabCombiner} object or \code{combinedTable}
#'
#' @param useID option to annotate identity-matched strings as "IDENTITY"
#'
#' @param minScore  numeric minimum allowable score (between 0 & 1) for
#'                  metabolomics feature pair alignments
#'
#' @param maxRankX  integer maximum allowable rank for X dataset features.
#'
#' @param maxRankY  integer maximum allowable rank for Y dataset features.
#'
#' @param method Conflict detection method. If equal to "score" (default),
#'              assigns a conflict subgroup if score of lower-ranking FPA is
#'              within some tolerance of higher-ranking FPA. If set to "mzrt",
#'              assigns a conflicting subgroup if within a small m/z & rt
#'              distance of the top-ranked FPA.
#'
#' @param delta numeric score or mz/rt distances used to define subgroups. If
#'              method = "score", a value (between 0 & 1) score difference
#'              between a pair of conflicting FPAs. If method = "mzrt", a length
#'              4 numeric: (m/z, rt, m/z, rt) tolerances, the first pair for X
#'              dataset features and the second pair for Y dataset features.
#'
#' @param maxRTerr numeric maximum allowable error between model-projected
#'                 retention time (rtProj) and observed retention time (rty)
#'
#' @param remove  Logical. Option to keep or discard rows deemed removable.
#'
#' @param resolveConflicts logical option to computationally resolve conflicting
#' rows to a final set of 1-1 feature pair alignments
#'
#' @param rtOrder logical. If resolveConflicts set to TRUE, then this imposes
#' retention order consistency on rows deemed "RESOLVED" within subgroups.
#'
#' @param balanced  Logical. Optional processing of "balanced" groups, defined
#'                  as groups with an equal number of features from input
#'                  datasets where all features have a 1-1 match.
#'
#' @param brackets_ignore character. If useID = TRUE, bracketed identity strings
#' of the types in this argument will be ignored
#'
#' @details
#' \code{metabCombiner} initially reports all possible feature pairings in the
#' rows of the \code{combinedTable} report. Most of these are misalignments that
#' require removal. This function is used to automate this reduction
#' process by labeling rows as removable or conflicting, based on certain
#' conditions, and is performed after computing similarity scores.
#'
#' A label may take on one of four values:
#'
#'  a) "": No determination made
#'  b) "IDENTITY": an alignment with matching identity "idx & idy" strings
#'  c) "REMOVE": a row determined to be a misalignment
#'  d) "CONFLICT": competing alignments for one or multiple shared features
#'
#' The labeling rules are as follows:
#'
#' 1) Groups determined to be 'balanced': label rows with rankX > 1 & rankY > 1
#'    "REMOVE" irrespective of \code{delta} criteria
#' 2) Rows with a score < \code{minScore}: label "REMOVE"
#' 3) Rows with rankX > \code{maxRankX} and/or rankY > \code{maxRankY}:
#'    label "REMOVE"
#' 4) Conflicting subgroup assignment as determined  by \code{method} &
#'    \code{delta} arguments. Conflicting alignments following outside
#'    \code{delta} thresholds: labeled "REMOVE". Otherwise, they are assigned
#'    a "CONFLICT" label and subgroup number.
#' 5) If \code{useID} argument set to TRUE, rows with matching idx & idy strings
#'    are labeled "IDENTITY". These rows are not changed to "REMOVE" or
#'    "CONFLICT" irrespective of subsequent criteria.
#'
#' @return  updated \code{combinedTable} or \code{metabCombiner} object. The
#' table will have three new columns:
#'
#' \item{labels}{characterization of feature alignments as described}
#' \item{subgroup}{conflicting subgroup number of feature alignments}
#' \item{alt}{alternate subgroup for rows in multiple feature pair conflicts}
#'
#' @examples
#'
#' #required steps prior to function use
#' data(plasma30)
#' data(plasma20)
#' p30 <- metabData(plasma30, samples = "CHEAR")
#' p20 <- metabData(plasma20, samples = "Red", rtmax = 17.25)
#' p.comb <- metabCombiner(xdata = p30, ydata = p20, binGap = 0.0075)
#' p.comb <- selectAnchors(p.comb, tolmz = 0.003, tolQ = 0.3, windy = 0.02)
#' p.comb <- fit_gam(p.comb, k = 20, iterFilter = 1)
#' p.comb <- calcScores(p.comb, A = 90, B = 14, C = 0.5)
#'
#' ##applies labels, but maintains all rows
#' p.comb <- labelRows(p.comb, maxRankX = 2, maxRankY = 2, maxRTerr = 0.5,
#'                     delta = 0.1, resolveConflicts = FALSE, remove = FALSE)
#'
#' ##automatically resolve conflicts and filter to 1-1 feature pairs
#' p.comb.2 <- labelRows(p.comb, resolveConflicts = FALSE, remove = FALSE)
#'
#' #this is identical to the previous command
#' p.comb.2 <- reduceTable(p.comb)
#'
#' p.comb = labelRows(p.comb, method = "mzrt", delta = c(0.005, 0.5, 0.005,0.3))
#'
#' ##this function may be applied to combinedTable inputs as well
#' cTable = cbind.data.frame(combinedTable(p.comb), featdata(p.comb))
#'
#' lTable = labelRows(cTable, maxRankX = 3, maxRankY = 2, minScore = 0.5,
#'          method = "score", maxRTerr = 0.5, delta = 0.2)
#'
#' @export
labelRows <- function(object, useID = FALSE, minScore = 0.5, maxRankX = 3,
                    maxRankY = 3, method = c("score", "mzrt"), delta = 0.1,
                    maxRTerr = 10, resolveConflicts = FALSE, rtOrder = TRUE,
                    remove = FALSE, balanced = TRUE,
                    brackets_ignore = c("(", "[", "{"))
{
    if(isCombinedTable(object) == 0) cTable <- object
    else if(isMetabCombiner(object) == 0) cTable <- combinedTable(object)
    else{
        combinerCheck(isMetabCombiner(object), "metabCombiner", "warning")
        combinerCheck(isCombinedTable(object), "combinedTable", "warning")
        stop("input object is not a valid metabCombiner or combinedTable")
    }
    method <- match.arg(method)
    check_lblrows_pars(maxRankX, maxRankY, minScore, maxRTerr, balanced,
                       method, delta, resolveConflicts, rtOrder)
    if(all(cTable[["score"]] == 1))
        stop("must call calcScores before using this function")
    fvs <- prepare_fields_values(cTable, useID, brackets_ignore)
    fields <- fvs[["fields"]]
    values <- fvs[["values"]]
    method <- as.integer(ifelse(method == "score", 1, 2))  #score: 1,  mzrt: 2
    rterr <- abs(fields[["rty"]] - fields[["rtProj"]])
    fields <- get_labels(fields, minScore, delta, maxRankX, maxRankY, method,
                        rterr, maxRTerr, balanced)
    if(resolveConflicts)  fields <- resolveRows(fields, rtOrder)
    cTable <- formLabeledTable(fields, values, remove, fvs[["group0"]])

    if(methods::is(object, "metabCombiner")){
        fdata <- featdata(object)
        fdata <- fdata[match(cTable[["rowID"]], fdata[["rowID"]]),]
        object <- update_mc(object, combinedTable = cTable, featdata = fdata)
    }
    else
        object <- cTable
    return(object)
}


#'@rdname labelRows
#'
#'@export
reduceTable <- function(object, useID = FALSE, maxRankX = 2, maxRankY = 2,
                        minScore = 0.5, delta = 0.1, maxRTerr = 10,
                        rtOrder = TRUE, brackets_ignore = c("(", "[", "{"))
{
    object <- labelRows(object, useID = useID, maxRankX = maxRankX,
                        maxRankY = maxRankY, minScore = minScore,
                        method = "score", delta = delta, maxRTerr = maxRTerr,
                        balanced = TRUE, resolveConflicts = TRUE,
                        rtOrder = rtOrder, remove = TRUE,
                        brackets_ignore = brackets_ignore)

    return(object)
}




