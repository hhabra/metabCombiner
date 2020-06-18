###Function for Evaluating A, B, C score Parameters

## piecewise function
hinge <- function(value, thresh){
    ifelse(value > thresh, value, 0)
}


#' @title Find Mismatch Identity Rows
#'
#' @description
#' Finds all rows containing misalignments for a given feature
#'
#' @param cTable  data frame. Abridged Combiner report table.
#'
#' @param id  character. Compound identity to be searched for misalignments
#'
#' @return vector of integer rows
#'
#' @noRd
mismatchfind <- function(cTable, id)
{
    mm = which(cTable[["idx"]] == id & cTable[["idy"]] != id |
               cTable[["idx"]] != id & cTable[["idy"]] == id)

    return(mm)
}

#' @title Calculate Mismatch Score
#'
#' @description Assign a mismatch value for a given metabolite, equal to the
#' highest mismatching combination for a given identified compound. If no
#' incorrect combinations exist, this is equal to 0.
#'
#' @param cTable data frame. Abridged metabCombiner report table.
#'
#' @param mismatches  integer rows containing misalignments
#'
#' @noRd
mismatchScore <- function(cTable, mismatches)
{
    ifelse(length(mismatches) > 0, max(cTable$score[mismatches]), 0)
}


#' @title Weight Parameter Objective Function
#'
#' @description
#' This function evaluates the A, B, C weight parameters in terms of score
#' separability of matching versus mismatching compound alignments. Higher
#' objective function value imply a superior weight parameter selection.
#'
#' @param cTable  data frame. Abridged \code{metabCombiner} report table.
#'
#' @param identities  data frame containing all evaluated identities
#'
#' @param A   Numeric weight for penalizing m/z differences.
#'
#' @param B   Numeric weight for penalizing differences between fitted &
#' observed retention times
#'
#' @param C Numeric weight for penalizing differences in Q (abundance quantiles)
#'
#' @param minScore numeric. Minimum score to count towards objective function value.
#'
#' @param mzdiff numeric differences between feature m/z values
#'
#' @param rtdiff Differences between model-projected retention time value &
#'               observed retention time
#'
#' @param qdiff  Difference between feature quantile Q values.
#'
#' @param rtrange  range of dataset Y retention times
#'
#' @param adductdiff Numeric divisors of computed score when non-empty adduct
#'                   labels do not match
#'
#' @param penalty  positive numeric penalty wherever S(i,j) > S(i,i), i =/= j
#'
#' @param matches  integer row indices of identity matches
#'
#' @param mismatches  list of integer mismatching identity rows for each identity
#'
#' @details
#' First, the similarity scores between all grouped features are calculated as
#' described in \code{scorePairs}
#'
#' Then, the objective value for a similarity S is evaluated as:
#'
#' \deqn{OBJ(S) = \sum h(S(i,i)) - h(S(i, j)) - p(S(i,i) > S(i,j))}
#'
#' -S(i,i) represents the similarity between correct identity alignments \cr
#' -S(i,j), represents the maximum similarity of i to grouped feature j,
#'         i =/= j (the highest-scoring misalignment) \cr
#' -h(x) = x if x > \code{minScore}, 0 otherwise \cr
#' -p(COND) = 0 if the condition is true, and a \code{penalty} value otherwise.
#'
#' This is summed over all labeled compound identities (e.g. idx = idy) shared
#' between input datasets.
#'
#' @return
#' A numeric value quantifying total separability of compound match similarity
#' scores from mismatch scores, given A,B,C values
objective <- function(cTable, identities, A, B, C, minScore, mzdiff, rtdiff,
                      qdiff, rtrange, adductdiff, penalty, matches, mismatches)
{
    cTable[["score"]] = scorePairs(A = A,
                                   B = B,
                                   C = C,
                                   mzdiff = mzdiff,
                                   rtdiff = rtdiff,
                                   qdiff = qdiff,
                                   rtrange = rtrange,
                                   adductdiff = adductdiff)

    identities[["match"]] = cTable[["score"]][matches]

    identities[["mismatch"]] = as.numeric(sapply(mismatches, function(mm)
                                                 mismatchScore(cTable, mm)))

    identities[["penalty"]]= ifelse(identities[["match"]] > identities[["mismatch"]],
                                     0, penalty)

    identities[["score"]] = hinge(identities[["match"]], minScore) -
                            hinge(identities[["mismatch"]], minScore) -
                            identities[["penalty"]]

    value = sum(identities[["score"]])

    return(value)
}

#' @title Evaluate Similarity Score Parameters
#'
#' @description
#' This function provides a method for guiding selection of suitable values for
#' A, B, & C weight arguments in the \code{\link{calcScores}} method, based on
#' the similarity scores of shared identified compounds. Datasets must have at
#' least one identity in common (i.e. idx = idy, case-insensitive), and
#' preferably more than 10.
#'
#' @param object metabCombiner object
#'
#' @param A  Numeric weights for penalizing m/z differences.
#'
#' @param B  Numeric weights for penalizing differences between fitted
#' & observed retention times
#'
#' @param C Numeric weights for penalizing differences in Q (abundance quantiles)
#'
#' @param fit Character. Choice of fitted rt model, "gam" or "loess."
#'
#' @param usePPM logical. Option to use relative ppm (as opposed to absolute) m/z
#'               differences in score computations.
#'
#' @param useAdduct logical. Option to penalize mismatches in (non-empty) adduct
#'                  column labels.
#'
#' @param adduct numeric. If useAdduct is TRUE, divides mismatching and non-empty
#'               adduct column labels by this value.
#'
#' @param minScore numeric. Minimum score to count towards objective function value.
#'
#' @param penalty  numeric. Subtractive mismatch penalty.
#'
#' @param groups integer. Vector of feature groups to score. If set to NULL
#'               (default), will compute scores for all feature groups.
#'
#' @param brackets_ignore bracketed identity character strings of these types
#' will be ignored according to this argument
#'
#' @return A data frame with the following columns:
#' \item{A}{m/z weight values}
#' \item{B}{rt weight values}
#' \item{C}{Q weight values}
#' \item{score}{objective function evaluation of (A,B,C) weights}
#'
#' @details
#' An objective function, based on the accurate and inaccurate alignments of shared
#' pre-identified compounds. For more details, see: \code{\link{objective}}.
#'
#' @note
#' In contrast to \code{\link{calcScores}} function, A, B, & C take numeric
#' vectors as input, as opposed to constants. The total number of rows in the
#' output will be equal to the products of the lengths of these input vectors
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
#' p.combined = selectAnchors(p.combined, tolMZ = 0.003, tolQ = 0.3, windY = 0.02)
#' p.combined = fit_gam(p.combined, k = seq(12,20,2), iterFilter = 1)
#'
#' ##example 1
#' scores = evaluateParams(p.combined, minScore = 0.7, penalty = 5)
#'
#' ##example 2: changing A, B, C parameters
#' scores = evaluateParams(p.combined, A = seq(70,150,5), B = seq(7, 20, 0.5),
#'                         C = c(0.2,0.3,0.5))
#'
#' ##example 3: using PPM mass deviation (note change to A argument)
#' scores = evaluateParams(p.combined, usePPM = TRUE, A = seq(0.01,0.1,0.01))
#'
#' ##example 4: limiting to groups 1-2000
#' scores = evaluateParams(p.combined, minScore = 0.5, groups = 1:2000)
#'
#' ##After evaluation
#' p.combined = calcScores(p.combined, A = 90, B = 14, C = 0.5)
#' }
#'
#' @seealso \code{\link{calcScores}}, \code{\link{objective}}
#'
#' @export
evaluateParams <- function(object, A = seq(60,150,by = 10), B = seq(6,15),
                       C = seq(0.1,0.5,by = 0.1) , fit = c("gam", "loess"),
                       usePPM = FALSE, useAdduct = FALSE, adduct = 1,
                       minScore = 0.5, penalty = 5, groups = NULL,
                       brackets_ignore = c("(", "[", "{"))
{
    code = isMetabCombiner(object)

    if(code)
        stop(combinerError(code, "metabCombiner"))

    cTable = combinedTable(object)[,1:15]

    if(any(is.na(A)) | any(is.na(B)) | any(is.na(C)))
        stop("At least one missing value in at least one of arguments A, B, C")

    if((class(A) != "numeric" & class(A) != "integer")|
       (class(A) != "numeric" & class(A) != "integer")|
       (class(A) != "numeric" & class(A) != "integer"))
        stop("arguments 'A', 'B', 'C' must be numeric vectors")

    if(any(A < 0) | any(B < 0) | any(C < 0))
        stop("arguments A, B, C must consist only of non-negative values")

    if(!is.numeric(minScore) | minScore < 0 | minScore > 1)
        stop("argument 'minScore' must be a positive constant between 0 & 1")

    if(!is.numeric(penalty) | penalty < 0)
        stop("argument 'penalty' must be a positive numeric constant")

    rtrange = max(cTable[["rty"]]) - min(cTable[["rty"]])

    #handling groups argument
    if(is.null(groups))
        groups = 1:max(cTable[["group"]])

    if(any(!groups %in% cTable[["group"]]))
        stop("Invalid argument 'groups'- at least one group value not detected")

    #rows whose scores will be computed according to groups
    rows = which(cTable[["group"]] %in% groups)
    cTable = cTable[rows,]

    cTable[["label"]] = compare_strings(cTable[["idx"]], cTable[["idy"]],
                                        "IDENTITY", "", brackets_ignore)

    if(!any(cTable[["label"]] == "IDENTITY"))
        stop("Must have at least one shared pre-labeled identity")

    cTable[["idx"]] = paste(tolower(cTable[["idx"]]),
                            cTable[["group"]],
                            sep = "_")

    cTable[["idy"]] = paste(tolower(cTable[["idy"]]),
                            cTable[["group"]],
                            sep = "_")

    identities = dplyr::filter(cTable, .data$label == "IDENTITY") %>%
                 dplyr::select(.data$idx) %>%
                 dplyr::mutate(`match` = 0,
                               `mismatch` = 0,
                               `penalty` = 0,
                               `score` = 0)

    names(identities)[1] = "id"

    counts = table(identities[["id"]])

    #look for duplicate identity names within a group
    if(any(counts > 1)){
        duplicates = which(counts > 1)
        dupliNames = names(counts)[duplicates]
        dupliGroups = unique(sapply(strsplit(dupliNames, "_"), "[[", 2))
        dupliGroups = paste(dupliGroups, collapse = ",")

        warning("At least one identity appears multiple times in a group.",
                "See groups: ", dupliGroups)

        identities = dplyr::filter(identities, !(id %in% duplicates))
    }

    ##reducing table to rows with identities of interest
    cTable = dplyr::filter(cTable,
                           .data$idx %in% identities[["id"]] |
                           .data$idy %in% identities[["id"]])

    fit = match.arg(fit)
    model = getModel(object, fit = fit)

    if(is.null(model))
        stop("object missing model of type ", fit)

    cTable[["rtProj"]] = predict(model, newdata = cTable)

    massdiff = mzdiff(mzx = cTable[["mzx"]],
                      mzy = cTable[["mzy"]],
                      usePPM = usePPM)

    rtdiff = abs(cTable[["rty"]] - cTable[["rtProj"]])
    qdiff = abs(cTable[["Qx"]] - cTable[["Qy"]])

    if(useAdduct == TRUE){
        if(!is.numeric(adduct) | adduct < 1){
            warning("value of argument 'adduct' must be a numeric value greater
                    than or equal to 1")
            adduct = 1
        }

        adductdiff = compare_adduct_strings(cTable[["adductx"]],
                                            cTable[["adducty"]],
                                            1, adduct, brackets_ignore)
    }
    else
        adductdiff = 1

    scores = data.frame(A = rep(A, each = length(B)*length(C)),
                        B = rep(B, each = length(C)),
                        C = C,
                        score = 0)

    matches = which(cTable[["label"]] == "IDENTITY")
    mismatches = lapply(identities[["id"]], function(id) mismatchfind(cTable,id))

    scores[["score"]] = mapply(function(A,B,C)
                               objective(cTable = cTable,
                                         identities = identities,
                                         A = A,
                                         B = B,
                                         C = C,
                                         minScore = minScore,
                                         mzdiff = massdiff,
                                         rtdiff = rtdiff,
                                         qdiff = qdiff,
                                         rtrange = rtrange,
                                         adductdiff = adductdiff,
                                         matches = matches,
                                         mismatches = mismatches,
                                         penalty = penalty),
                                scores[["A"]],
                                scores[["B"]],
                                scores[["C"]])

    scores = dplyr::arrange(scores, desc(.data$score))

    return(scores)
}
