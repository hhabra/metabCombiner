

## Determine Matching Identity Strings
#'
#' @param s1  One of two character vectors to be compared
#'
#' @param s2  One of two character vectors to be compared
#'
#' @param value The value to give for matching s1 and s2 strings
#'
#' @return  character vector. Indices of matching strings are equal to value,
#'          with remaining indices left as empty characters.
#'
#' @noRd
##
comp_id_strings <- function(s1,s2, value){
    s1 = ifelse(is.na(s1), "", s1)
    s2 = ifelse(is.na(s2), "", s2)

    ifelse(!(s1 == "" | s2 == "") & tolower(s1) == tolower(s2), value, "")
}


##
#' @title Label and Remove Combiner Report Rows
#'
#' @description
#' Updates the \code{combinerTable} with labeled determinations of which rows are
#' removable, as well as conflicting and identity-labeled rows. Optionally, these
#' removable rows may be discarded, resulting in a reduced \code{combinerTable}.
#'
#' @param object  Either a \code{metabCombiner} object or \code{combinerTable}.
#'
#' @param maxRankX  Integer. Maximum allowable rank for X dataset features.
#'                  Must be greater than 1.
#'
#' @param maxRankY  Integer. Maximum allowable rank for Y dataset features.
#'                  Must be greater than 1.
#'
#' @param minScore  Numeric. Minimum allowable score (between 0 & 1) for
#'                  metabolomics feature matches.
#'
#' @param balanced  Logical. Optional processing of "balanced" groups (defined as
#'                  groups with equal numbers of features from X & Y datasets).
#'                  If no contradicting top-matches,
#'
#' @param conflict numeric vector. Non-negative Tolerance values for determining
#'                 if pairs of features within a single dataset are considered
#'                 "conflicting. Must be a length 4 vector. The first and second
#'                 values give the m/z and rt tolerances, respectively, for
#'                 dataset X; third and fourth value give the m/z tolerance and rt
#'                 tolerances, respectively, for dataset Y.
#'
#' @param remove  Logical. Option to keep or discard rows deemed removable (labeled
#'                "REMOVE" in Combiner Report)
#'
#' @details
#' Combiner initially reports all possible feature alignments in the rows of
#' \code{combinerTable} report. Most of these alignments are inaccurate and
#' require inspection and removal. This function is used to automate most of the
#' reduction process by labeling rows as removable or conflicting, based on
#' certain conditions, and is performed after computing similarity scores.
#'
#' The labeling rules are as follows:
#' 1) Rows with matching idx & idy strings are labeled "IDENTITY". These rows are
#' not labeled "REMOVE", irrespective of subsequent criteria.
#' 2) Rows with a score < \code{minScore}: label "REMOVE"
#' 3) Rows with rankX > \code{maxRankX} or rankY > \code{maxRankY}: label "REMOVE"
#' 4) If row has score > \code{minScore} and 1 < rankX \eqn{\leq} \code{maxRankX}:
#'    or 1 < rankY \eqn{\leq}\code{maxRankY}:
#'
#' @return  combinerTable with label column for removable rows
#'
#' @export
##
reduceTable <- function(object, maxRankX = 3, maxRankY = 3, minScore = 0.3,
                        balanced = TRUE, conflict = c(0.003,0.3, 0.003,0.3),
                        remove = FALSE)
{
    if(isCombinerTable(object) == 0)
        cTable = object

    else{
        code = isMetabCombiner(object)

        if(code)
            stop(combinerError(code, "metabCombiner"))

        cTable = combinerTable(object)
    }

    if(!class(maxRankX) %in% c("numeric", "integer") |
       !class(maxRankY) %in% c("numeric", "integer"))
        stop("arguments maxRankX & maxRankY must be numeric")

    maxRankX = as.integer(maxRankX[1])
    maxRankY = as.integer(maxRankY[1])

    if(maxRankX < 1 | maxRankY < 1)
        stop("arguments maxRankX & maxRankY must be greater than 1")

    if(class(minScore) != "numeric" | minScore > 1 | minScore < 0)
        stop("argument 'minScore' must be a numeric value between 0 and 1")

    if(!is.logical(balanced))
        stop("expected a logical for argument 'balanced'")

    if(!is.numeric(conflict) | length(conflict) != 4)
        stop("argument 'conflict' must be a length four numeric vector")

    if(any(conflict < 0))
        stop("tolerance values in parameter 'conflict' must be non-negative")

    cTable = cTable[with(cTable, order(group, desc(score))), ]

    values = cTable[,-c(1:15)]
    fields = cTable[,1:15]
    rm(cTable)

    #bugfix for existing columns named "labels"
    if(any(names(values) == "labels")){
        remove = which(names(values) == "labels")
        values = values[,-remove]
    }

    fields[["labels"]] = comp_id_strings(fields[["idx"]], fields[["idy"]],
                                         "IDENTITY")

    fields[["labels"]] = .Call("findRemovables",
                               labels = fields[["labels"]],
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
                               PACKAGE = "Combiner")

    #option to eliminate rows labeled as removables
    if(remove == TRUE){
        keepRows = which(fields[["labels"]] != "REMOVE")
        fields = fields[keepRows,]
        values = values[keepRows,]
    }

    cTable = data.frame(fields, values, stringsAsFactors = FALSE,
                        check.names = FALSE)

    if(class(object) == "metabCombiner")
        object@combinerTable = cTable

    else
        object = cTable

    return(object)
}


