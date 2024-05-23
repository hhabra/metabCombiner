##
#' @title Calculate m/z Differences
#'
#' @description
#' Helper function for \code{calcScores} & \code{evaluateParams}. Calculates
#' the absolute or relative (ppm) m/z differences between grouped features.
#'
#' @param mzx numeric vector. m/z values from dataset X.
#'
#' @param mzy  numeric vector. m/z values from dataset Y.
#'
#' @param usePPM logical. Option to calculate relative differences (ppm).
#'
#' @return A vector of positive numeric relative or absolute m/z differences.
#'
#' @noRd
##
massdiff <- function(mzx, mzy, usePPM)
{
    diffs <- mzx - mzy

    if(usePPM){
        mins <- ifelse(diffs < 0, mzx, mzy)
        diffs <- 1e6*abs(diffs)/ mins
    }

    return(diffs)
}


#' Parameterize m/z Shift Model
#'
#' @description Return a list of parameters for computing a m/z shift model
#' between a pair of LC-MS metabolomics data sets.
#'
#' @param mz Numeric. Limits the m/z distance between feature pairs for modeling
#'
#' @param rt Numeric. Limits the RT distance between feature pairs for modeling
#' (note: this is proportion of total retention time)
#'
#' @param Q Numeric. Limits the Q distance between feature pairs for modeling
#'
#' @param k integer k value controlling the dimension of the basis spline fit
#'
#' @details Correcting for systematic m/z shifts improves the scores for
#' feature pair alignments, yielding more accurate match hypotheses. This
#' function yields a parameter list for GAM spline fitting of points that
#' meet criteria for m/z, RT (rty vs rtProj), and Q tolerances. The number
#' of knots for the GAM fit is controlled by hyperparameter k.
#'
#' @return list of parameter values
#'
#' @export
#'
#' @examples
#'
#' mzfit <- mzfitParam(mz = 0.003, rt = 0.05, Q = 0.2, k = 20)
#'
mzfitParam <- function(mz = 0.003, rt = 0.03, Q = 0.3, k = 10)
{
    return(list(mz = mz, rt = rt, Q = Q, k = k))
}



#' Compute m/z Shift Model
#'
#' @description Generates a GAM model for correcting systematic m/z shifts
#' observed between a pair of LC-MS data sets.
#'
#' @param object metabCombiner object
#'
#' @param fit List of m/z shift parameter values; see ?mzfitParam
#'
#' @param plot Logical. Option to plot the m/z shift model.
#'
#' @param ... other arguments to be passed to plot
#'
#' @details Correcting for systematic m/z shifts improves the scores for
#' feature pair alignments, yielding more accurate match hypotheses. This
#' function generates a basis spline curve, modeling the m/z shift (mzy - mzx)
#' as a function of m/z (mzx). Selected points are ordered feature pairs that
#' meet criteria for m/z, RT (rty vs rtProj), and Q tolerances set in the
#' \code{mzfit} list argument (see: \code{\link{mzfitParam}}). If
#'
#' Setting the \code{plot} option to TRUE generates an image of the curve fit
#' through the selected points and is a useful method for determining if
#' m/z mapping is appropriate for the analysis and tuning certain parameters.
#'
#' This function is called within \code{\link{calcScores}}, which can help
#' improve the pairwise scores
#'
#' @return model object of class gam or 0 (if no model selected)
#'
#' @export
#'
#' @examples
#'
#' \donttest{
#' data(plasma30)
#' data(plasma20)
#'
#' p30 <- metabData(plasma30, samples = "CHEAR")
#' p20 <- metabData(plasma20, samples = "Red", rtmax = 17.25)
#' p.comb <- metabCombiner(xdata = p30, ydata = p20, binGap = 0.0075)
#'
#' p.comb <- selectAnchors(p.comb, tolmz = 0.003, tolQ = 0.3, windy = 0.02)
#' p.comb <- fit_gam(p.comb, k = 20, iterFilter = 1, family = "gaussian")
#'
#' mzmodel <- mzFit(p.comb, plot = TRUE,
#'                  fit = mzfitParam(mz = 0.003, rt = 0.03, Q = 0.3, k = 20))
#'
#'}
#'
mzFit <- function(object, fit = mzfitParam(), plot = TRUE, ...)
{
    if(is(object, "metabCombiner")){
        cTable <- combinedTable(object)
        if(is.null(getModel(object))) stop("No RT fit detected!")
        cTable[["rtProj"]] <- predict(getModel(object), newdata = cTable)
    }
    else
        cTable <- object
    if(all(cTable[["rtProj"]] == 0)){
        warning("retention time mapping must precede m/z fitting")
        return(mgcv::gam(0 ~ s(mzx)), data = rTable)
    }

    rTable <- cTable[combinerNames()]
    rtrange <- max(rTable$rty, na.rm = TRUE) - min(rTable$rty, na.rm = TRUE)
    rTable <- dplyr::filter(rTable, abs(.data$mzy - .data$mzx) < fit[["mz"]],
                        abs(.data$rty - .data$rtProj) < fit[["rt"]] * rtrange,
                        abs(.data$Qy - .data$Qx) < fit[["Q"]])
    model <- mgcv::gam((mzy - mzx) ~ s(mzx, k = fit[["k"]], bs = "tp"),
                       data = rTable)
    if(plot){
        plot(rTable[["mzx"]], rTable[["mzy"]] - rTable[["mzx"]],
             main = "m/z Shift Model", xlab = "m/z", ylab = "m/z shift",
             pch = 18, ...)
        lines(rTable[["mzx"]], predict(model, rTable), lwd = 3, col = "red")
    }

    return(model)
}




#' @title Calculate Pairwise Alignment Scores
#'
#' @description Helper function for \code{\link{calcScores}} &
#' \code{\link{evaluateParams}}. Calculates a pairwise similarity score between
#' grouped features using differences in m/z, rt, and Q.
#'
#' @param A Numeric weight for penalizing m/z differences.
#'
#' @param B Numeric weight for penalizing differences between fitted & observed
#' retention times.
#'
#' @param C Numeric weight for differences in Q (abundance quantiles).
#'
#' @param mzdiff Numeric differences between feature m/z values
#'
#' @param rtdiff Differences between model-projected retention time value &
#'              observed retention time
#'
#' @param qdiff Difference between feature quantile Q values
#'
#' @param rtrange  Range of dataset Y retention times
#'
#' @param adductdiff Numeric divisors of computed score when non-empty adduct
#'                  labels do not match
#'
#' @details
#' The score between two grouped features x & y is calculated as:
#'
#' \deqn{S = -exp(-A |mzx - mzy| - B |rty - rtproj|/rtrange - C |Qx - Qy|)}
#'
#' where \code{mzx} & \code{Qx} correspond to the m/z and abundance quantile
#' values of feature x; \code{mzy}, \code{rty}, and \code{Qy} correspond to the
#' m/z, retention time, and quantile values of feature y; \code{rtproj}
#' is the model-projected retention time of feature x onto the Y dataset
#' chromatogram and \code{rtrange} is the retention time range of the Y dataset
#' chromatogram. \code{A}, \code{B}, \code{C} are non-negative constant weight
#' parameters for penalizing m/z, rt, and Q  differences. Values between 0 (no
#' confidence alignment) and 1 (high confidence alignment).
#'
#' @return  Numeric similarity score between 0 & 1
scorePairs <- function(A, B, C, mzdiff, rtdiff, qdiff, rtrange, adductdiff)
{
    nlogScore <- A * abs(mzdiff) + B * abs(rtdiff) / rtrange + C * abs(qdiff)
    score <- base::exp(-nlogScore) / adductdiff
    return(score)
}

#' Rank and Arrange Alignments by Score
#'
#' @description  Rank feature pair alignments within each group in the
#' \code{combinedTable} by their score. This is limited to rows specified
#' within the main \code{\link{calcScores}} function.
#'
#' @param cTable \code{combinedTable} of metabCombiner object
#'
#' @param rows integer rows to apply the operations in \code{combinedTable}
#'
#' @return \code{combinedTable} with updated \code{rankX} & \code{rankY}
#' columns arranged by group and (descending) score order
#'
#' @noRd
calculateRanks <- function(cTable, rows)
{
    ##bug fix for duplicated column names
    cTable_copy <- cTable[, c("mzx", "mzy", "rtx", "rty", "score", "rankX", "rankY")]

    cTable_copy[rows,] <- cTable_copy[rows,] %>%
        dplyr::group_by(.data$mzx, .data$rtx) %>%
        dplyr::mutate(rankX = dplyr::dense_rank(desc(.data$score))) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(.data$mzy, .data$rty) %>%
        dplyr::mutate(rankY = dplyr::dense_rank(desc(.data$score))) %>%
        dplyr::ungroup()

    cTable[c("rankX", "rankY")] = cTable_copy[c("rankX", "rankY")]
    cTable <- cTable[with(cTable,order(`group`, desc(`score`))),]

    return(cTable)
}


##
#' @title Compute Feature Similarity Scores
#'
#' @description Calculates a pairwise similarity (between 0 & 1) between all
#' grouped features in \code{metabCombiner} object. The similarity score
#' calculation is described in \code{\link{scorePairs}}.
#'
#' @param object  metabCombiner object.
#'
#' @param A Numeric weight for penalizing m/z differences.
#'
#' @param B Numeric weight for penalizing differences between fitted & observed
#' retention times
#'
#' @param C Numeric weight for differences in Q (abundance quantiles).
#'
#' @param fit Character. Choice of fitted rt model, "gam" or "loess."
#'
#' @param mzshift Logical. If TRUE, shifts the m/z values (mzx) before scoring.
#'
#' @param mzfit List of parameters for shifting m/z values; see ?mzfitParam
#'
#' @param usePPM logical. Option to use relative (as opposed to absolute) m/z
#'              differences in score computations.
#'
#' @param useAdduct logical. Option to penalize mismatches in (non-empty,
#'                  non-bracketed) adduct column annotations.
#'
#' @param adduct numeric. If useAdduct is TRUE, divides score of mismatched,
#'               non-empty and non-bracked adduct column labels by this value.
#'
#' @param groups integer. Vector of feature groups to score. If set to NULL
#'              (default), will compute scores for all feature groups.
#'
#' @param brackets_ignore If useAdduct = TRUE, bracketed adduct character
#' strings of these types will be ignored according to this argument
#'
#' @return \code{metabCombiner} object with updated \code{combinedTable}.
#' rtProj column will contain fitted retention times determined from previously
#' computed model; score will contain computed pairwise similarity scores of
#' feature pairs; rankX & rankY are the integer ranks of scores for x & y
#' features in descending order.
#'
#' @details
#' This function updates the \code{rtProj}, \code{score}, \code{rankX}, and
#' \code{rankY} columns in the \code{combinedTable} report. First, using the
#' RT mapping model computed in the previous steps, \code{rtx} values are
#' projected onto \code{rty}. Then similarity scores are calculated based on
#' m/z, rt (rtProj vs rty), and Q differences, with multiplicative weight
#' penalties \code{A}, \code{B}, and \code{C}.
#'
#' If the datasets contain representative set of shared identities (idx = idy),
#' \code{\link{evaluateParams}} provides some guidance on appropriate \code{A},
#' \code{B}, and \code{C} values to use. In testing, the best values for
#' \code{A} should lie between 50 and 120, according to mass accuracy; if using
#' ppm (\code{usePPM} = TRUE), the suggested range is between 0.01 and 0.05.
#' \code{B} should be between 5 and 15 depending on fitting accuracy (higher
#' if datasets processed under roughly identical conditions) ; \code{C} should
#' vary between 0 and 1, depending on sample similarity. See examples below.
#'
#' Some input datasets exhibit systematic m/z shifts
#'
#'  If using adduct information (\code{useAdduct} = TRUE), the score is divided
#'  by the numeric \code{adduct} argument if non-empty and non-bracketed adduct
#'  values do not match. Be sure that adduct annotations are accurate before
#'  using this functionality.
#'
#' @seealso
#' \code{\link{evaluateParams}}, \code{\link{scorePairs}}
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
#' p.comb <- selectAnchors(p.comb, tolmz = 0.003, tolQ = 0.3, windy = 0.02)
#' p.comb <- fit_gam(p.comb, k = 20, iterFilter = 1, family = "gaussian")
#'
#' #example: moderate m/z deviation, accurate rt fit, high sample similarity
#' p.comb <- calcScores(p.comb, A = 90, B = 14, C = 0.8, useAdduct = FALSE,
#'          groups = NULL, fit = "gam", usePPM = FALSE)
#' cTable <- combinedTable(p.comb)  #to view results
#'
#' \donttest{
#' #example 2: high m/z deviation, moderate rt fit, low sample similarity
#' p.comb <- calcScores(p.comb, A = 50, B = 8, C = 0.2)
#'
#' #example 3: low m/z deviation, poor rt fit, moderate sample similarity
#' p.comb <- calcScores(p.comb, A = 120, B = 5, C = 0.5)
#'
#' #example 4: using ppm for mass deviation; note different A value
#' p.comb <- calcScores(p.comb, A = 0.05, B = 14, C = 0.5, usePPM = TRUE)
#'
#' #example 5: limiting to specific m/z groups 1-1000
#' p.comb <- calcScores(p.comb, A = 90, B = 14, C = 0.5, groups = seq(1,1000))
#'
#' #example 6: using adduct information
#' p.comb <- calcScores(p.comb, A = 90, B = 14, C = 0.5, useAdduct = TRUE,
#'                      adduct = 1.25)
#' }
#' @export
##
calcScores <- function(object, A = 75, B = 10, C = 0.25, groups = NULL,
                        fit = c("gam", "loess"), mzshift = FALSE,
                        mzfit = mzfitParam(), useAdduct = FALSE, adduct = 1.25,
                        usePPM = FALSE, brackets_ignore = c("(", "[", "{"))
{
    combinerCheck(isMetabCombiner(object), "metabCombiner")
    cTable <- combinedTable(object)
    fdata <- featData(object)
    fit <- match.arg(fit)
    model <- getModel(object, fit = fit)
    if(is.null(groups)) groups <- setdiff(unique(cTable[["group"]]),0)
    check_score_pars(cTable, A, B, C, fit, model, groups, adduct = adduct)

    if(length(A) > 1 | length(B) > 1 | length(C) > 1){
        warning("only the first element used for arguments 'A', 'B', 'C'")
        A <- A[1];   B <- B[1];   C <- C[1]
    }

    rtrange <- max(cTable$rty, na.rm = TRUE) - min(cTable$rty, na.rm = TRUE)
    rows <-  which(cTable[["group"]] %in% setdiff(groups, 0))
    cTable$rtProj[rows] <- stats::predict(model, newdata = cTable[rows,])
    if(useAdduct){
        if(!is.numeric(adduct) | adduct < 1)
            stop("'adduct' argument must be a (>=1) numeric value")
        adductdiff <- compare_strings(cTable$adductx[rows],cTable$adducty[rows],
                                    1, adduct, brackets_ignore, type = "mm")
    }
    else
        adductdiff <- 1

    if(isTRUE(mzshift))
        mzS <- stats::predict(mzFit(object, mzfit, FALSE), newdata = cTable)
    else
        mzS <- rep(0, length(cTable$mzx))

    mzdiff <-  massdiff(cTable$mzx[rows] + mzS[rows], cTable$mzy[rows], usePPM)
    rtdiff <- abs(cTable$rty[rows] - cTable$rtProj[rows])
    qdiff <- abs(cTable$Qx[rows] - cTable$Qy[rows])

    cTable$score[rows] <- scorePairs(A = A, B = B, C = C, mzdiff = mzdiff,
                                    rtdiff = rtdiff, qdiff = qdiff,
                                    rtrange = rtrange, adductdiff = adductdiff)

    cTable[c("rtProj", "score")] <- round(cTable[c("rtProj", "score")], 4)
    cTable <- calculateRanks(cTable, rows)
    fdata <- fdata[match(cTable[["rowID"]], fdata[["rowID"]]),]

    object <- update_mc(object, combinedTable = cTable, featData = fdata,
                        coefficients = list(`A` = A, `B` = B, `C` = C))
    return(object)
}


