##
#' @title Calculate m/z Differences
#'
#' @description
#' Helper function for \code{calcScores} & \code{evaluateParams}. Calculates the
#' absolute or relative (ppm) m/z differences between grouped features.
#'
#' @param mzx  numeric vector. m/z values from dataset X.
#'
#' @param mzy  numeric vector. m/z values from dataset Y.
#'
#' @param usePPM logical. Option to calculate relative differences (ppm).
#'
#' @return A vector of positive numeric relative or absolute m/z differences.
#'
#' @noRd
##
mzdiff <- function(mzx, mzy, usePPM){
    diffs = abs(mzx - mzy)

    if(usePPM)
        diffs = 1e6*diffs/ sapply(1:length(mzx), function(r) min(mzx[r],mzy[r]))

    return(diffs)
}

##
#' @title Determine if two adduct strings match.
#'
#' @description
#' Helper function for \code{calcScores} & \code {evaluateParams}.Determine if two
#' strings (feature adduct labels) match. If both strings are non-empty and not
#' equivalent, return adduct penalty; if matching strings or at least one of the
#' two is empty,return 1. Penalty values will divide the feature similarity score.
#'
#' @param s1  One of two character vectors to be compared
#'
#' @param s2  One of two character vectors to be compared
#'
#' @param adduct Numeric (value >= 1) multiplicative adduct mismatch penalty.
#'
#' @return Numeric vector equal to 1 or user-supplied mismatch penalty.
#'
#' @noRd
##
comp_adduct_strings <- function(s1,s2, adduct){
    s1 = ifelse(is.na(s1), "", s1)
    s2 = ifelse(is.na(s2), "", s2)

    ifelse(!(s1 == "" | s2 == "") & tolower(s1) != tolower(s2),
         adduct, 1)
}

##
#' @title Calculate the pairwise score between grouped features
#'
#' @description Helper function for \code{calcScores()} & \code{evalParams()}.
#' Calculates the pairwise similarity score between grouped features acquired
#' in complementary datasets.
#'
#' @param A Numeric weight for penalizing m/z differences.
#'
#' @param B Numeric weight for penalizing differences between fitted & observed
#' retention times
#'
#' @param C Numeric weight for differences in Q (abundance quantiles).
#'
#' @param mzdiff Numeric differences between feature m/z values
#'
#' @param rtdiff Differences between model-projected retention time value &
#'               observed retention time
#'
#' @param qdiff Difference between feature quantile Q values
#'
#' @param rtrange   Range of dataset Y retention times
#'
#' @param adductdiff Numeric divisors of computed score when non-empty adduct
#'                   labels do not match
#'
#' @details
#' Score between two grouped features x & y calculated as:
#'
#' \deqn{S(x,y,A,B,C) = -exp(-A |mzx-mzy|- B |rty-rtproj|/rtrange - C |Qx-Qy|)}
#'
#' where mzx & Qx correspond to the m/z and abundance quantile values of feature
#' x, mzy, rty, and Qy correspond to the m/z, retention time, and abundance
#' quantile values of feature y,  rtproj is the model-projected retention of
#' feature x onto the retention times of the dataset containing y, and rtrange is
#' the range of retention times of y dataset features. A, B, C are non-negative
#' constant weight parameters for penalizing m/z, rt, and Q  differences. Values
#' vary between 0 (no confidence combination) and 1 (high confidence combination).
#'
#'
#' @return  Pairwise score between two grouped features(between 0 & 1) evaluating
#'          the likelihood of a compound match.
#'
#'
scorePairs <- function(A, B, C, mzdiff, rtdiff, qdiff, rtrange, adductdiff){

    nlogScore = A * abs(mzdiff) + B * abs(rtdiff) / rtrange + C * abs(qdiff)

    score = base::exp(-nlogScore) / adductdiff

    return(score)
}


##
#' @title Compute Feature Similarity Scores
#'
#' @description Calculates a pairwise similarity (between 0 & 1) between all
#' grouped features in \code{metabCombiner} object.
#'
#' @param object  metabCombiner object.
#'
#' @param A Numeric weight for penalizing m/z differences.
#'
#' @param B Numeric weight for penalizing differences between fitted & observed
#' retention times
#'
#' @param C Numeric weight for penalizing differences in Q (abundance quantiles).
#'
#' @param fit Character. Choice of fitted rt model, "gam" or "loess."
#'
#' @param usePPM logical. Option to use relative (as opposed to absolute) m/z
#'               differences in score computations.
#'
#' @param useAdduct logical. Option to penalize mismatches in (non-empty) adduct
#'                  column labels.
#'
#' @param adduct numeric. If useAdduct is TRUE, divides mismatching and non-empty
#'               adduct column labels by this value.
#'
#' @param groups integer. Vector of feature groups to score. If set to NULL
#'              (default), will compute scores for all feature groups.
#'
#' @return metabCombiner object with updated combinerTable. rtProj column will
#' contain fitted retention times determined from previously computed model;
#' score will contain the computed pairwise similarity scores of features from
#' datasets x & y; rankX & rankY are the integer ranks of scores for x & y
#' features in descending order.
#'
#' @export
##
calcScores <- function(object, A, B, C, fit = c("gam", "loess"), usePPM = FALSE,
                       useAdduct = FALSE, adduct = 2, groups = NULL)
{
    code = isMetabCombiner(object)

    if(code)
        stop(combinerError(code, "metabCombiner"))

    if(class(A) != "numeric" | class(B) != "numeric" | class(C) != "numeric")
        stop("arguments 'A', 'B', 'C' must be numeric constants")

    if(length(A) > 1 | length(B) > 1 | length(C) > 1){
        warning("Expected constants for arguments 'A', 'B', 'C' and
                only the first element will be used.")

        A = A[1];   B = B[1];   C = C[1]
    }

    if(A < 0 | B < 0 | C < 0)
        stop("arguments 'A', 'B', 'C' must be non-negative")

    cTable = combinerTable(object)

    fit = match.arg(fit)
    model = getModel(object, fit = fit)

    if(is.null(model))
        stop(paste("object missing model of type ", fit, sep =""))

    rtrange = max(cTable[["rty"]]) - min(cTable[["rty"]])

    if(useAdduct == FALSE)
        adduct = 1

    #handling groups argument
    if(is.null(groups))
        groups = 1:max(cTable[["group"]])

    if(any(!groups %in% cTable[["group"]]))
        stop("Invalid argument 'groups'- at least one group value not detected")

    #rows whose scores are computed according to groups
    rows = which(cTable[["group"]] %in% groups)

    cTable$rtProj[rows] = predict(model, newdata = cTable[rows,])

    massdiff = mzdiff(mzx = cTable$mzx[rows],
                      mzy = cTable$mzy[rows],
                      usePPM = usePPM)

    rtdiff = abs(cTable$rty[rows] - cTable$rtProj[rows])
    qdiff = abs(cTable$Qx[rows] - cTable$Qy[rows])

    adductdiff = comp_adduct_strings(cTable$adductx[rows],
                                     cTable$adducty[rows],
                                     adduct = adduct)

    cTable$score[rows] = scorePairs(A = A,
                                    B = B,
                                    C = C,
                                    mzdiff = massdiff,
                                    rtdiff = rtdiff,
                                    qdiff = qdiff,
                                    rtrange = rtrange,
                                    adductdiff = adductdiff
                                    )

    cTable[c("rtProj", "score")] = round(cTable[c("rtProj", "score")], 4)

    ##bug fix for duplicated column names
    cT = cTable[ , c("mzx", "mzy", "rtx", "rty", "score", "rankX", "rankY")]

    ##calculate score ranking of features
    cT[rows,] =   cT[rows,] %>%
                  group_by(mzx, rtx) %>%
                  mutate(rankX = dense_rank(desc(score))) %>%
                  ungroup() %>%
                  group_by(mzy, rty) %>%
                  mutate(rankY= dense_rank(desc(score))) %>%
                  ungroup()

    cTable[c("rankX", "rankY")] = cT[c("rankX", "rankY")]

    object@combinerTable = cTable[with(cTable, order(group, desc(score))), ]
    object@coefficients = list(A = A, B = B, C = C)

    return(object)
}


