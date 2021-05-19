
#' Reform Report with New Labels
#'
#' @description Helper function for labelRows(). After determining row
#' annotations, stitches together the metadata, the row annotations, and
#' the sample + extra value data.
#'
#' @param fields data.frame abridged combinedTable with metadata fields.
#'
#' @param values data.frame combinedTable with sample (+ extra) columns
#'
#' @noRd
formLabeledTable <- function(fields, values, remove)
{
    #option to eliminate rows labeled as removables
    if(remove == TRUE){
        keepRows <- which(fields[["labels"]] != "REMOVE")
        fields <- fields[keepRows,]
        values <- values[keepRows,]
    }

    #when there are existing columns named "labels", "subgroup", "alt"
    labnames <- c("labels", "subgroup", "alt")
    N <- 1

    while(any(names(values) %in% labnames)){
        labnames <- paste(c("labels", "subgroup", "alt"), N, sep = ".")
        names(fields)[seq(16,18)] <- labnames
        N <- N+1
    }

    cTable <- data.frame(fields, values, stringsAsFactors = FALSE,
                        check.names = FALSE)

    return(cTable)
}


#' @title Annotate and Remove Report Rows
#'
#' @description
#' This is a method for annotating identity-matched, removable, &
#' conflicting feature pair alignment (FPA) rows in the \code{combinedTable}
#' report. Simple thresholds for score, rank, retention time error and delta
#' score can computationally reduce the set of possible FPAs to the most likely
#' compound matches. FPAs falling within some small measure (in score or mz/rt)
#' of the top-ranked row are organized into subgroups to facilitate inspection;
#' setting delta to 0 automatically reduces to 1-1 matches.
#'
#' \code{reduceTable} behaves identically to labelRows, but with delta set to 0
#' & remove set to TRUE, automatically limiting to 1 - 1 feature matches
#' constrained by rank and score threshold parameters. Rank thresholds defaults
#' are also stricter with reduceTable.
#'
#' @param object Either a \code{metabCombiner} object or \code{combinedTable}.
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
#' @param balanced  Logical. Optional processing of "balanced" groups, defined
#'                  as groups with an equal number of features from input
#'                  datasets where all features have a 1-1 match.
#'
#' @param remove  Logical. Option to keep or discard rows deemed removable.
#'
#' @param brackets_ignore character. Bracketed identity strings of the types
#' in this argument will be ignored
#'
#' @details
#' \code{metabCombiner} initially reports all possible FPAs in the rows of the
#' \code{combinedTable} report. Most of these are misalignments that
#'  require removal. This function is used to automate most of the reduction
#'  process by labeling rows as removable or conflicting, based on certain
#'  conditions, and is performed after computing similarity scores.
#'
#'  A label may take on one of four values:
#'
#'  a) "": No determination made
#'  b) "IDENTITY": an alignment with matching identity "idx & idy" strings
#'  c) "REMOVE": a row determined to be a misalignment
#'  d) "CONFLICT": competing alignments for one or multiple shared features
#'
#' The labeling rules are as follows:
#' 1) Rows with matching idx & idy strings are labeled "IDENTITY". These rows
#'    are not labeled "REMOVE", irrespective of subsequent criteria.
#' 2) Groups determined to be 'balanced': label rows with rankX > 1 & rankY > 1
#'    "REMOVE" irrespective of \code{delta} criteria
#' 3) Rows with a score < \code{minScore}: label "REMOVE"
#' 4) Rows with rankX > \code{maxRankX} and/or rankY > \code{maxRankY}:
#'    label "REMOVE"
#' 5) Conflicting subgroup assignment as determined  by \code{method} &
#'    \code{delta} arguments. Conflicting alignments following outside
#'    \code{delta} thresholds: labeled "REMOVE". Otherwise, they are assigned
#'    a "CONFLICT" label and subgroup number.
#'
#' @return  updated \code{combinedTable} or \code{metabCombiner} object. The
#' table will have three new columns:
#'
#' \item{labels}{characterization of feature alignments as described}
#' \item{subgroup}{conflicting subgroup number of feature alignments}
#' \item{alt}{alternate subgroup for rows in multiple feature pair conflicts}
#'
#' @examples
#' data(plasma30)
#' data(plasma20)
#'
#' p30 <- metabData(plasma30, samples = "CHEAR")
#' p20 <- metabData(plasma20, samples = "Red", rtmax = 17.25)
#' p.comb = metabCombiner(xdata = p30, ydata = p20, binGap = 0.0075)
#' p.comb = selectAnchors(p.comb, tolmz = 0.003, tolQ = 0.3, windy = 0.02)
#' p.comb = fit_gam(p.comb, k = 20, iterFilter = 1)
#' p.comb = calcScores(p.comb, A = 90, B = 14, C = 0.5)
#'
#' ###merge combinedTable and featdata
#' cTable = cbind.data.frame(combinedTable(p.comb), featdata(p.comb))
#'
#' ##example using score-based conflict detection method
#' lTable = labelRows(cTable, maxRankX = 3, maxRankY = 2, minScore = 0.5,
#'          method = "score", maxRTerr = 0.5, delta = 0.2)
#'
#' ##example using mzrt-based conflict detection method
#' lTable = labelRows(cTable, method = "mzrt", maxRankX = 3, maxRankY = 2,
#'                      delta = c(0.005, 0.5, 0.005,0.3), maxRTerr = 0.5)
#'
#' @export
labelRows <- function(object, minScore = 0.5, maxRankX = 3, maxRankY = 3,
                    method = c("score", "mzrt"), delta = 0.1, maxRTerr = 10,
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
                       method, delta)
    cTable <- cTable[with(cTable, order(`group`, desc(`score`))), ]
    values <- cTable[,-seq(1,15)]
    fields <- cTable[,seq(1,15)]
    fields[["labels"]] <- compare_strings(fields[["idx"]], fields[["idy"]],
                                        "IDENTITY", "", brackets_ignore)
    fields[["subgroup"]] <- fields[["alt"]]  <- integer(nrow(fields))
    method <- as.integer(ifelse(method == "score", 1, 2))  #score: 1,  mzrt: 2
    rterr <- abs(fields[["rty"]] - fields[["rtProj"]])

    fields[["labels"]] <- .Call("labelRows",
                            labels = fields$labels, subgroup = fields$subgroup,
                            alt = fields$alt, mzx = fields$mzx,
                            mzy = fields$mzy, rtx = fields$rtx,
                            rty = fields$rty, score = fields$score,
                            rankX = fields$rankX, rankY = fields$rankY,
                            group = fields$group, balanced = balanced,
                            delta = delta, minScore = minScore,
                            maxRankX = as.integer(maxRankX[1]),
                            maxRankY = as.integer(maxRankY[1]),
                            method = method, maxRTerr = as.numeric(maxRTerr),
                            rterr = rterr, PACKAGE = "metabCombiner")

    cTable <- formLabeledTable(fields, values, remove)

    if(methods::is(object, "metabCombiner")){
        fdata <- featdata(object)[row.names(cTable),]
        object <- update_mc(object, combinedTable = cTable, featdata = fdata)
    }
    else
        object <- cTable

    return(object)
}


#'@rdname labelRows
#'
#'@export
reduceTable <- function(object, maxRankX = 2, maxRankY = 2, minScore = 0.5,
                        maxRTerr = 10, brackets_ignore = c("(", "[", "{"))
{
    object <- labelRows(object, maxRankX = maxRankX, maxRankY = maxRankY,
                        minScore = minScore, method = "score", delta = 0,
                        maxRTerr = maxRTerr, balanced = TRUE,
                        brackets_ignore = brackets_ignore)

    return(object)
}




