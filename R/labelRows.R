
#' @title Label and Remove metabCombiner Report Rows
#'
#' @description
#' Automated method for determining identity-matching, removable, and conflicting
#' alignments in \code{combinedTable}. Optionally, removable rows may be
#' eliminated, resulting in a reduced \code{combinedTable}. Remaining alignments
#' that require further consideration are labeled "CONFLICT" and organized into
#' subgroups if they fall within some small measure (in score or mz/rt).
#'
#' @param object  Either a \code{metabCombiner} object or \code{combinedTable}.
#'
#' @param maxRankX  Integer. Maximum allowable rank for X dataset features.
#'
#' @param maxRankY  Integer. Maximum allowable rank for Y dataset features.
#'
#' @param minScore  Numeric. Minimum allowable score (between 0 & 1) for
#'                  metabolomics feature alignments.
#'
#' @param method  Conflict detection method. If equal to "score" (default method),
#'                assigns a conflict subgroup if score of lower-ranking alignment
#'                is within some tolerance of higher-ranking alignment. If set to
#'                "mzrt", assigns a conflicting subgroup if within a small m/z &
#'                rt distance of the top-ranked alignment.
#'
#' @param conflict Numeric. If method = "score", a constant (between 0 & 1)
#'                 score difference between a pair of alignments.
#'                 If method = "mzrt", a length 4 numeric: (m/z, rt, m/z, rt)
#'                 tolerances, the first two for X dataset features and the
#'                 second two for Y dataset features.
#'
#' @param balanced  Logical. Optional processing of "balanced" groups (defined as
#'                  groups with equal numbers of features from X & Y datasets with
#'                  no conflicting top-matches).
#'
#' @param remove  Logical. Option to keep or discard rows deemed removable (labeled
#'                "REMOVE" in \code{combinedTable})
#'
#' @param brackets_ignore character. Bracketed identity strings of the types
#' in this argument will be ignored
#'
#' @details
#' \code{metabCombiner} initially reports all possible feature alignments in the
#'  rows of \code{combinedTable} report. Most of these are misalignments that
#'  require inspection and removal. This function is used to automate most of the
#'  reduction process by labeling rows as removable or conflicting, based on
#'  certain conditions, and is performed after computing similarity scores.
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
#'    "REMOVE" irrespective of \code{conflict} criteria
#' 3) Rows with a score < \code{minScore}: label "REMOVE"
#' 4) Rows with rankX > \code{maxRankX} and/or rankY > \code{maxRankY}:
#'    label "REMOVE"
#' 5) Conflicting subgroup assignment as determined  by \code{method} &
#'    \code{conflict} arguments. Conflicting alignments following outside
#'    \code{conflict} thresholds: labeled "REMOVE". Otherwise,
#'
#' @return  updated \code{combinedTable} or \code{metabCombiner} object. The
#' tabel will have three new columns:
#'
#' \item{labels}{characterization of feature alignments as described}
#' \item{subgroup}{conflicting subgroup number of feature alignments}
#' \item{alt}{alternative subgroup for rows conflicting with multiple top-matches}
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
#' p.combined = selectAnchors(p.combined, tolMZ = 0.003, tolQ = 0.3, windY = 0.02)
#' p.combined = fit_gam(p.combined, k = seq(12,20,2), iterFilter = 1)
#' p.combined = calcScores(p.combined, A = 90, B = 14, C = 0.5)
#' cTable = combinedTable(p.combined)
#'
#' ##example using score-based conflict detection method
#' lTable = labelRows(cTable, method = "score", conflict = 0.2)
#'
#' ##example using mzrt-based conflict detection method
#' lTable = labelRows(cTable, method = "mzrt", conflict = c(0.005, 1, 0.005,0.5))
#'
#' ##example changing thresholds for rank and minimum score
#' lTable = labelRows(cTable, minScore = 0.5, maxRankX = 2, maxRankY = 2,
#'                    method = "score", conflict = 0.1)
#' }
#' @export
labelRows <- function(object, maxRankX = 3, maxRankY = 3, minScore = 0.3,
                      conflict, method = c("score", "mzrt"), balanced = TRUE,
                      remove = FALSE, brackets_ignore = c("(", "[", "{"))
{
    if(iscombinedTable(object) == 0)
        cTable = object

    else{
        code = isMetabCombiner(object)

        if(code)
            stop(combinerError(code, "metabCombiner"))

        cTable = combinedTable(object)
    }

    if(!class(maxRankX) %in% c("numeric", "integer") |
      !class(maxRankY) %in% c("numeric", "integer"))
        stop("arguments maxRankX & maxRankY must be integers")

    maxRankX = as.integer(maxRankX[1])
    maxRankY = as.integer(maxRankY[1])

    if(maxRankX < 1 | maxRankY < 1)
        stop("arguments maxRankX & maxRankY must be greater than 1")

    if(class(minScore) != "numeric" | minScore > 1 | minScore < 0)
        stop("argument 'minScore' must be a numeric value between 0 and 1")

    if(!is.logical(balanced))
        stop("expected a logical for argument 'balanced'")

    method = match.arg(method)

    if(method == "score"){
        if(!is.numeric(conflict) | conflict > 1 | conflict < 0)
            stop("argument 'conflict' must be a numeric value between 0 & 1")

        if(length(conflict) > 1)
            stop("with method = score, constant value expected for conflict")
    }

    else if(method == "mzrt"){
        if(length(conflict) != 4)
            stop("with method = mzrt, length 4 vector expected for conflict)")

        if(any(conflict < 0))
            stop("values in 'conflict' argument must be non-negative")
    }

    cTable = cTable[with(cTable, order(`group`, desc(`score`))), ]

    values = cTable[,-c(1:15)]
    fields = cTable[,1:15]
    rm(cTable)

    fields[["labels"]] = compare_strings(fields[["idx"]], fields[["idy"]],
                                       "IDENTITY", "", brackets_ignore)

    fields[["subgroup"]] = integer(nrow(fields))
    fields[["alt"]] = integer(nrow(fields))

    #score: 1,  mzrt: 2
    method = as.integer(ifelse(method == "score", 1, 2))

    fields[["labels"]] = .Call("labelRows",
                                labels = fields[["labels"]],
                                subgroup = fields[["subgroup"]],
                                alt = fields[["alt"]],
                                mzx = fields[["mzx"]],
                                mzy = fields[["mzy"]],
                                rtx = fields[["rtx"]],
                                rty = fields[["rty"]],
                                score = fields[["score"]],
                                rankX = fields[["rankX"]],
                                rankY = fields[["rankY"]],
                                group = fields[["group"]],
                                balanced = balanced,
                                conflict = conflict,
                                minScore = minScore,
                                maxRankX = maxRankX,
                                maxRankY = maxRankY,
                                method = method,
                                PACKAGE = "metabCombiner")

      #option to eliminate rows labeled as removables
      if(remove == TRUE){
          keepRows = which(fields[["labels"]] != "REMOVE")
          fields = fields[keepRows,]
          values = values[keepRows,]
      }

    #when there are existing columns named "labels", "subgroup", "alt"
    labnames = c("labels", "subgroup", "alt")
    N = 1

    while(any(names(values) %in% labnames)){
        labnames = paste(c("labels", "subgroup", "alt"), N, sep = ".")
        names(fields)[16:18] = labnames
        N = N+1
    }

    cTable = data.frame(fields, values, stringsAsFactors = FALSE,
                        check.names = FALSE)

    if(class(object) == "metabCombiner")
        object@combinedTable = cTable

    else
        object = cTable

    return(object)
}


