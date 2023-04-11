###Function for Evaluating A, B, C score Parameters

## piecewise function
hinge <- function(value, thresh){
    ifelse(value > thresh, value, 0)
}

#' Form Evaluation Table for Matched IDs
#'
#' @param cTable data.frame. Modified combinerTable containing only metadata
#' and labels for matching strings.
#'
#' @return data.frame with the following columns:
#' \item{id}{names of shared identified known features}
#' \item{match}{scores of correctly matched alignments of known features}
#' \item{mismatch}{maximum mismatch score for respective id}
#' \item{penalty}{penalty value if mismatch score > match score; 0 otherwise}
#' \item{score}{score = match - mismatch - penalty}
#'
#' @noRd
form_idtable <- function(cTable){
    idtable <- dplyr::filter(cTable, .data$label == "IDENTITY") %>%
        dplyr::select(.data$idx) %>%
        dplyr::mutate(`match`= 0, `mismatch`= 0,`penalty`= 0,`score`= 0) %>%
        dplyr::rename("id" = "idx")

    #look for duplicate identity names within a group
    counts <- table(idtable[["id"]])

    if(any(counts > 1)){
        duplicates <- which(counts > 1)
        dupliNames <- names(counts)[duplicates]
        dupliGroups <- paste(unique(vapply(strsplit(dupliNames, "_"), "[[", 2,
                                    FUN.VALUE = character(1))), collapse = ",")
        warning("At least one identity appears multiple times in a group.",
                "See groups: ", dupliGroups)
        idtable <- dplyr::filter(idtable, !(id %in% duplicates))
    }

    if(nrow(idtable) == 0)
        stop("no valid shared identities for evaluation")

    return(idtable)
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
    mm <- which(cTable[["idx"]] == id & cTable[["idy"]] != id |
                cTable[["idx"]] != id & cTable[["idy"]] == id)

    return(mm)
}

#' @title Calculate Mismatch Score
#'
#' @description Assign a mismatch value for a given metabolite, equal to the
#' highest mismatching combination for a given identified compound. If no
#' incorrect combinations exist, this is equal to 0.
#'
#' @param cTable data.frame. Abridged metabCombiner report table.
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
#' @param idtable  data frame containing all evaluated identities
#'
#' @param A   Numeric weight for penalizing m/z differences.
#'
#' @param B   Numeric weight for penalizing differences between fitted &
#' observed retention times
#'
#' @param C Numeric weight for differences in Q (abundance quantiles).
#'
#' @param minScore numeric. Minimum score to count towards objective value.
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
#' @param mismatches list of integer identity row mismatches for each identity
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
#' -p(COND) = 0 if the condition is true, and a \code{penalty} value otherwise
#'
#' This is summed over all labeled compound identities (e.g. idx = idy) shared
#' between input datasets.
#'
#' @return
#' A numeric value quantifying total separability of compound match similarity
#' scores from mismatch scores, given A,B,C values
objective <- function(cTable, idtable, A, B, C, minScore, mzdiff, rtdiff,
                    qdiff, rtrange, adductdiff, penalty, matches, mismatches)
{
    cTable[["score"]] <- scorePairs(A = A,
                                    B = B,
                                    C = C,
                                    mzdiff = mzdiff,
                                    rtdiff = rtdiff,
                                    qdiff = qdiff,
                                    rtrange = rtrange,
                                    adductdiff = adductdiff)

    idtable[["match"]] <- cTable[["score"]][matches]

    idtable[["mismatch"]] <- as.numeric(vapply(mismatches, function(mm)
                                    mismatchScore(cTable, mm), numeric(1)))

    idtable[["penalty"]] <- ifelse(idtable[["match"]] > idtable[["mismatch"]],
                                    0, penalty)

    idtable[["score"]] <- hinge(idtable[["match"]], minScore) -
                            hinge(idtable[["mismatch"]], minScore) -
                            idtable[["penalty"]]

    value <- sum(idtable[["score"]])

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
#' @param C Numeric weight for differences in Q (abundance quantiles).
#'
#' @param fit Character. Choice of fitted rt model, "gam" or "loess."
#'
#' @param usePPM logical. Option to use relative parts per million (ppm) as
#'              opposed to absolute) m/z differences in score computations.
#'
#' @param minScore numeric minimum score to count towards objective function
#' calculation for known matching features (idx = idy) and mismatches.
#'
#' @param penalty numeric. Subtractive mismatch penalty.
#'
#' @param groups integer. Vector of feature groups to score. If set to NULL
#'              (default), will compute scores for all feature groups.
#'
#' @param brackets_ignore bracketed identity and adduct character strings of
#' these types will be ignored according to this argument
#'
#' @return A data frame with the following columns:
#' \item{A}{m/z weight values}
#' \item{B}{rt weight values}
#' \item{C}{Q weight values}
#' \item{totalScore}{objective function evaluation of (A,B,C) weights}
#'
#' @details
#' This uses an objective function, based on the accurate and inaccurate
#' alignments of  shared pre-identified compounds. For more details, see:
#' \code{\link{objective}}.
#'
#' @note
#' In contrast to \code{\link{calcScores}} function, A, B, & C take numeric
#' vectors as input, as opposed to constants. The total number of rows in the
#' output will be equal to the products of the lengths of these input vectors
#'
#' @examples
#'
#' data(plasma30)
#' data(plasma20)
#'
#' p30 <- metabData(plasma30, samples = "CHEAR")
#' p20 <- metabData(plasma20, samples = "Red", rtmax = 17.25)
#' p.comb <- metabCombiner(xdata = p30, ydata = p20, binGap = 0.0075)
#'
#' p.comb <- selectAnchors(p.comb, windx = 0.03, windy = 0.02)
#' p.comb <- fit_gam(p.comb, k = 20, iterFilter = 2)
#'
#' #example 1
#' scores <- evaluateParams(p.comb, A = seq(60,100,10), B = seq(10,15), C = 0.5,
#'     minScore = 0.7, penalty = 10)
#'
#' ##example 2: limiting to groups 1-2000
#' scores <- evaluateParams(p.comb, minScore = 0.5, groups = seq(1,2000))
#'
#' @seealso \code{\link{calcScores}}, \code{\link{objective}}
#'
#' @export
evaluateParams <- function(object, A = seq(60,150,by = 10), B = seq(6,15),
                        C = seq(0.1,0.5,by = 0.1) , fit = c("gam", "loess"),
                        usePPM = FALSE, minScore = 0.5, penalty = 5,
                        groups = NULL,brackets_ignore = c("(", "[", "{"))
{
    combinerCheck(isMetabCombiner(object), "metabCombiner")
    cTable <- combinedTable(object)[combinerNames()]
    fit <- match.arg(fit)
    model <- getModel(object, fit = fit)
    if(is.null(groups))  groups <- setdiff(unique(cTable[["group"]]),0)
    check_score_pars(cTable, A, B, C, model, fit, groups, minScore, penalty)
    rtrange <- max(cTable$rty, na.rm = TRUE) - min(cTable$rty, na.rm = TRUE)
    rows <- which(cTable[["group"]] %in% groups)
    cTable <- cTable[rows,]
    cTable[["label"]] <- compare_strings(cTable[["idx"]], cTable[["idy"]],
                                        "IDENTITY", "", brackets_ignore)
    if(!any(cTable[["label"]] == "IDENTITY"))
        stop("must have at least one shared labeled identities (idx = idy)")

    cTable[["idx"]] <- paste(tolower(cTable[["idx"]]), cTable[["group"]],
                            sep = "_")
    cTable[["idy"]] <- paste(tolower(cTable[["idy"]]),cTable[["group"]],
                            sep = "_")
    idtable <- form_idtable(cTable)
    cTable <- dplyr::filter(cTable, .data$idx %in% idtable[["id"]] |
                            .data$idy %in% idtable[["id"]])

    cTable[["rtProj"]] <- predict(model, newdata = cTable)
    massdiff <- mzdiff(cTable[["mzx"]], cTable[["mzy"]], usePPM)
    rtdiff <- abs(cTable[["rty"]] - cTable[["rtProj"]])
    qdiff <- abs(cTable[["Qx"]] - cTable[["Qy"]])
    scores <- data.frame(A = rep(A, each = length(B)*length(C)),
                        B = rep(B, each = length(C)), C = C, score = 0)

    matches <- which(cTable[["label"]] == "IDENTITY")
    mismatches <- lapply(idtable[["id"]], function(id) mismatchfind(cTable,id))

    scores[["totalScore"]] <- mapply(function(A,B,C)
                    objective(cTable = cTable, idtable = idtable,
                                A = A, B = B, C = C, minScore = minScore,
                                mzdiff = massdiff, rtdiff = rtdiff,
                                qdiff = qdiff, rtrange = rtrange,
                                adductdiff = 1, matches = matches,
                                mismatches = mismatches, penalty = penalty),
                                scores[["A"]],  scores[["B"]], scores[["C"]])

    scores <- dplyr::arrange(scores, desc(.data$totalScore))
    return(scores)
}
