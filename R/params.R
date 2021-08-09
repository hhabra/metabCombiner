#' @title List selectAnchors Defaults
#'
#' @description
#' List of default parameters for anchor selection step of main package
#' workflow, which can be used as input for the wrapper functions. See
#' help(selectAnchors) or ?selectAnchors for more details.
#'
#' @param useID   Choice of using IDs for anchor selection; default: FALSE
#' @param tolmz   m/z tolerance for ordered pair features; default: 0.003
#' @param tolQ    Q tolerance for ordered pair features; default: 0.3
#' @param tolrtq  RT quantile tolerance for ordered pair features; default: 0.5
#' @param windx   X feature RT window parameter. Default: 0.03
#' @param windy   Y feature RT window parameter. Default: 0.03
#' @param brackets_ignore   bracket types for ignoring string comparisons
#'
#' @return list of selectAnchors parameters
#'
#' @examples
#' sa_param <- selectAnchorsParam(tolmz = 0.002, tolQ = 0.2, windy = 0.02)
#'
#' @seealso \code{\link{selectAnchors}}, \code{\link{metabCombine}}
#'
#' @export
selectAnchorsParam <- function(useID = FALSE, tolmz = 0.003, tolQ = 0.3,
                                tolrtq = 0.3, windx = 0.03, windy = 0.03,
                                brackets_ignore = c("(", "[", "{"))
{
    params <- list(useID = useID, brackets_ignore = brackets_ignore,
                   tolmz = tolmz, tolQ = tolQ, tolrtq = tolrtq,
                   windx = windx, windy = windy)

    return(params)
}

#' @title List fit_gam Defaults
#'
#' @description
#' List of default parameters for GAM fitting step of main package workflow,
#' which can be used as input for the wrapper functions. See help(fit_gam)
#' or ?fit_gam for more details.
#'
#' @param useID     choice of preserving identity-based anchors; default: FALSE
#' @param k           values for GAM basis dimension k
#' @param iterFilter  number of outlier filtering iterations; default: 2
#' @param outlier   outlier filtering method (either "MAD" (mean absolute
#'                  deviation) or "boxplot"); default: "MAD"
#' @param coef    outlier filtering coefficient; default: 2
#' @param prop  minimum proportion of fits in which a point can be a flagged
#'              outlier; default: 0.5
#' @param weights   optional supplied weights to individual points; default: 1
#' @param bs      choice of spline type ("bs" or "ps"); default: "bs"
#' @param family  choice of family ("scat" or "gaussian"); default: "scat"
#' @param m       basis and penalty order; default: c(3,2)
#' @param method  smoothing parameter estimation method; default: "REML"
#' @param optimizer  numerical optimization for GAM; default: "newton"
#' @param message   option to print progress message; default: TRUE
#'
#' @return list of fit_gam parameters
#'
#' @examples
#' fitParam <- fitgamParam(k = c(12,14,18,20), iterFilter = 1, bs = "ps",
#'                         family = "gaussian", method = "GCV.Cp")
#'
#' @seealso \code{\link{fit_gam}}, \code{\link{metabCombine}}
#'
#' @export
fitgamParam <- function(useID = FALSE, k = seq(10,20,2), iterFilter = 2,
                        outlier = "MAD", coef = 2, prop = 0.5, weights = 1,
                        bs = "bs", family = "scat", m = c(3,2), method = "REML",
                        optimizer = "newton", message = TRUE)
{
    params <- list(useID = useID, k = k, iterFilter = iterFilter,
                   outlier = outlier, coef = coef, prop = prop, bs = bs,
                   family = family, weights = weights, m = m, method = method,
                   optimizer = optimizer, message = message)
    return(params)
}



#' @title List fitLoess Defaults
#'
#' @description
#' List of default parameters for loess fitting step of main package workflow,
#' See help(fit_loess) or ?fit_loess for more details.
#'
#' @return list of fit_loess parameters:
#'
#' @param useID choice of preserving identity-based anchors; default: FALSE
#' @param spans values for span parameter which controls degree of smoothing
#' @param iterFilter number of outlier filtering iterations; default: 2
#' @param outlier outlier filtering method (either "MAD" or "boxplot");
#'      default: "MAD"
#' @param coef outlier filtering coefficient; default: 2
#' @param prop minimum proportion of fits where a point can be a flagged outlier;
#'       default: 0.5
#' @param weights optional supplied weights to individual points; default: 1
#' @param message option to print progress message; default: TRUE
#' @param control loess-specific control parameters; see: ?loess.control
#'
#' @seealso \code{\link{fit_loess}}, \code{\link{metabCombine}}
#'
#' @examples
#' fitParam <- fitloessParam(spans = c(0.2,0.25,0.3), outlier = "boxplot",
#'                          iterFilter = 3, coef = 1.5, message = FALSE,
#'                          control = loess.control(iterations = 4))
#'
#' @export
fitloessParam <- function(useID = FALSE, spans = seq(0.2, 0.3, by = 0.02),
                        outlier = "MAD", coef = 2, iterFilter = 2,
                        prop = 0.5, weights = 1, message = TRUE,
                        control = loess.control(surface = "direct",
                                                iterations = 10))
{
    params <- list(useID = useID, spans = spans, outlier = outlier, coef = coef,
                    iterFilter = iterFilter, prop = prop, weights = weights,
                    message = message, control = control)
    return(params)
}

#' @title List calcScores Defaults
#'
#' @description
#' List of default parameters for score calculation step of main package
#' workflow. See help(calcScores) or ?calcScores for details.
#'
#' @return list of calcScores parameters
#'
#' @param A m/z difference specific weight; default: 75
#' @param B RT prediction error specific weight; default: 10
#' @param C Q difference specific weight; default: 0.25
#' @param fit choice of fitted model ("gam" or "loess"); default: "gam"
#' @param groups choice of m/z groups to score
#' @param usePPM choice to use PPM for m/z differences; default: FALSE
#' @param useAdduct choice to use adduct strings in scoring; default: FALSE
#' @param adduct value divisor for mismatched adduct strings; default: 1.25
#' @param brackets_ignore bracket types for ignoring string comparisons
#'
#' @examples
#' cs_param <- calcScoresParam(A = 60, B = 15, C = 0.3)
#'
#' cs_param <- calcScoresParam(A = 0.1, B = 20, C = 0.2, usePPM = TRUE)
#'
#' @seealso \code{\link{calcScores}}, \code{\link{metabCombine}}
#'
#' @export
calcScoresParam <- function(A = 75, B = 10, C = 0.25, fit = "gam",
                            groups = NULL, usePPM = FALSE, useAdduct = FALSE,
                            adduct = 1.25, brackets_ignore = c("(", "[", "{"))
{
    params <- list(A = A, B = B, C = C, fit = fit, groups = groups,
                   usePPM = usePPM, useAdduct = useAdduct, adduct = adduct,
                   brackets_ignore = brackets_ignore)
    return(params)

}

#' @title List labelRows & reduceTable Defaults
#'
#' @description
#' List of default parameters for combinedTable row annotation and removal.
#' See help(labelRows) or ?labelRows for more details. reduceTableParam loads
#' parameters for the more automated \code{reduceTable} function
#'
#' @return list of labelRows parameters
#'
#' @param maxRankX  maximum rank allowable for X features
#' @param maxRankY  maximum rank allowable for Y features
#' @param minScore  minimum score threshold; default: 0.5
#' @param method thresholding method for subgroup detection ("score" or "mzrt");
#'      default: "score"
#' @param delta  score distance or mz/rt difference tolerances for subgrouping;
#'      default: 0.1
#' @param maxRTerr maximum allowable difference between predicted RT (rtProj) &
#'     observed RT (rty); default: 10 minutes
#' @param resolveConflicts logical. If TRUE, automatically resolves subgroups to
#'     1-1 feature pair alignments
#' @param rtOrder logical. If TRUE and resolveConflicts is TRUE, imposes
#'     retention order condition on paired alignments
#'
#' @param remove option to eliminate rows determined as removable;
#'     default: FALSE
#' @param balanced option to reduce balanced groups; default: TRUE

#' @param brackets_ignore bracket types for ignoring string comparisons
#'
#' @examples
#' lrParams <- labelRowsParam(maxRankX = 2, maxRankY = 2, delta = 0.1,
#'                              maxRTerr = 0.5)
#'
#' @seealso \code{\link{labelRows}}, \code{\link{metabCombine}},
#' \code{\link{batchCombine}}, \code{\link{reduceTable}}
#'
#' @export
labelRowsParam <- function(maxRankX = 3, maxRankY = 3, minScore = 0.5,
                            delta = 0.1, method = "score", maxRTerr = 10,
                            resolveConflicts = FALSE, rtOrder = TRUE,
                            remove = FALSE,  balanced = TRUE,
                            brackets_ignore = c("(", "[", "{"))
{
    params <- list(maxRankX = maxRankX, maxRankY = maxRankY,
                  minScore = minScore, method = method, delta = delta,
                  maxRTerr = maxRTerr, balanced = balanced, remove = remove,
                  resolveConflicts = resolveConflicts, rtOrder = rtOrder,
                  brackets_ignore = brackets_ignore)
    return(params)
}

#' @rdname labelRowsParam
#'
#' @export
reduceTableParam <- function(maxRankX = 2, maxRankY = 2, minScore = 0.5,
                            maxRTerr = 10, delta = 0.1, rtOrder = TRUE,
                            brackets_ignore = c("(", "[", "{"))
{
    params <- list(maxRankX = maxRankX, maxRankY = maxRankY,
                   minScore = minScore, method = "score", delta = delta,
                   maxRTerr = maxRTerr, balanced = TRUE, remove = TRUE,
                   resolveConflicts = TRUE, rtOrder = rtOrder,
                   brackets_ignore = brackets_ignore)
    return(params)
}

