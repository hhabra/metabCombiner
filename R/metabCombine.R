#' @title metabCombiner Wrapper Function
#'
#' @description \code{metabCombine} wraps the five main \code{metabCombiner}
#' workflow steps into a single wrapper function. Parameter list arguments
#' organize program parameters by constituent package functions.
#'
#' @param xdata metabData object. One of two datasets to be combined.
#'
#' @param ydata metabData object. One of two datasets to be combined.
#'
#' @param binGap numeric parameter used for grouping features by m/z.
#' See ?mzGroup for more details.
#'
#' @param xid character identifier of xdata. If xdata is a metabData, assigns a
#' new ID for this dataset; if xdata is a metabCombiner, must be assigned to one
#' of the existing dataset IDs. See details for more information.
#'
#' @param yid character identifier of ydata. If ydata is a metabData, assigns a
#' new ID for this dataset; if ydata is a metabCombiner, must be assigned to one
#' of the existing dataset IDs. See details for more information.
#'
#' @param means logical. Option to take average m/z, rt, and/or Q from
#' \code{metabComber}. May be a vector (length = 3), single value (TRUE/FALSE),
#' or a list with names "mz", "rt", "Q" as names.
#'
#' @param fitMethod RT spline-fitting method, either "gam" or "loess"
#'
#' @param union logical. Option to include non-matched features in final
#' \code{combinedTable} results
#'
#' @param impute logical. If TRUE, imputes the mean mz/rt/Q values for missing
#' features in \code{metabCombiner} object inputs before use in alignment (not
#' recommended for disparate data alignment); if FALSE, features with missing
#' information are dropped.
#'
#' @param anchorParam list of parameter values for selectAnchors() function
#'
#' @param fitParam list of parameter values for fit_gam() or fit_loess()
#'
#' @param scoreParam list of parameter values for calcScores()
#'
#' @param labelParam list of parameter values for labelRows()
#'
#' @param rtOrder logical. If set to TRUE, retention order consistency expected
#' when resolving conflicting alignments for \code{metabCombiner} object inputs.
#'
#' @details
#' The five main steps in \code{metabCombine} are 1) m/z grouping & combined
#' table construction, 2) selection of ordered pair RT anchors, 3) nonlinear
#' spline (Basis Spline GAM or LOESS) fitting to predict RTs, 4) score
#' calculation and feature pair alignment ranking, 5) combined table row
#' annotation and reduction. metabData arguments \code{xdata} & \code{ydata}
#' and m/z grouping \code{binGap} are required for step 1.
#'
#' Steps 2-5 are handled by \code{anchors}, \code{fit}, \code{scores}, &
#' \code{labels}, respectively, with lists containing the argument values for
#' each step expected for these arguments. \code{\link{selectAnchorsParam}},
#' \code{\link{fitgamParam}}, \code{\link{fitloessParam}},
#' \code{\link{calcScoresParam}}, & \code{\link{labelRowsParam}} load the
#' default program values of \code{\link{selectAnchors}}, \code{\link{fit_gam}},
#' \code{\link{fit_loess}}, \code{\link{calcScores}} & \code{\link{labelRows}},
#' respectively. These program arguments should be modified as necessary for
#' the datasets used for analysis.
#'
#' By default, the RT fitting method (\code{fitMethod}) is set to "gam", which
#' means the argument \code{fit} is a list of parameters for \code{fit_gam};
#' if the (\code{fitMethod}) argument is set to "loess", then the \code{fit}
#' argument expects a list of \code{fit_loess} parameters.
#'
#' @return a \code{metabCombiner} object following complete analysis
#'
#' @seealso \code{\link{selectAnchorsParam}}, \code{\link{fitgamParam}},
#' \code{\link{calcScoresParam}}, \code{\link{labelRowsParam}},
#' \code{\link{fitloessParam}}
#'
#' @examples
#' \donttest{
#' data("plasma20")
#' data("plasma30")
#'
#' p30 <- metabData(plasma30, samples = "CHEAR")
#' p20 <- metabData(plasma20, samples = "Red", rtmax = 17.25)
#'
#' #parameter lists:
#' saParam <- selectAnchorsParam(tolrtq = 0.2, windy = 0.02, tolmz = 0.002)
#' fitParam <- fitgamParam(k = seq(12,15), iterFilter = 1, outlier = "boxplot",
#'                         family = "gaussian", prop = 0.6, coef = 1.5)
#' scoreParam <- calcScoresParam(A = 75, B = 15, C = 0.3)
#' labelParam <- labelRowsParam(maxRankX = 2, maxRankY = 2, delta = 0.1)
#'
#' #metabCombine wrapper
#' p.combined <- metabCombine(xdata = p30, ydata = p20, binGap = 0.0075,
#'                            anchorParam = saParam, fitParam = fitParam,
#'                            scoreParam = scoreParam, labelParam = labelParam)
#'
#' ##to view results
#' p.combined.table <- combinedTable(p.combined)
#'
#' }
#'
#' @export
metabCombine <- function(xdata, ydata, binGap = 0.005, xid = NULL, yid = NULL,
                        means = list('mz' = FALSE, 'rt' = FALSE, 'Q' = FALSE),
                        fitMethod = "gam", rtOrder = TRUE, union = FALSE,
                        impute = FALSE, anchorParam = selectAnchorsParam(),
                        fitParam = fitgamParam(),scoreParam = calcScoresParam(),
                        labelParam = labelRowsParam())
{
    anchors <- anchorParam;    fit <- fitParam;
    scores <- scoreParam;      labels <- labelParam;
    object <- metabCombiner(xdata = xdata, ydata = ydata, binGap = binGap,
                            xid = xid, yid = yid, means = means,
                            rtOrder = rtOrder, impute = impute)
    object <- selectAnchors(object, useID = anchors$useID, tolQ = anchors$tolQ,
                            tolmz = anchors$tolmz, tolrtq = anchors$tolrtq,
                            windx = anchors$windx, windy = anchors$windy,
                            brackets_ignore = anchors$brackets_ignore)

    if(fitMethod == "gam")
        object <- fit_gam(object, useID = fit$useID, k = fit$k,
                            iterFilter = fit$iterFilter, outlier = fit$outlier,
                            coef = fit$coef, prop = fit$prop, bs = fit$bs,
                            family = fit$family, m = fit$m, method = fit$method,
                            rtx = fit$rtx, rty = fit$rty, weights = fit$weights,
                            optimizer = fit$optimizer, message = fit$message)
    else if(fitMethod == "loess")
        object <- fit_loess(object, useID = fit$useID, spans = fit$spans,
                            iterFilter = fit$iterFilter, outlier = fit$outlier,
                            coef = fit$coef, prop = fit$prop, rtx = fit$rtx,
                            rty = fit$rty, weights = fit$weights,
                            message = fit$message)
    else
        stop("'fitMethod' must be either 'gam' or 'loess'")

    object <- calcScores(object, A = scores$A, B = scores$B, C = scores$C,
                        fit = fitMethod, groups = scores$groups,
                        usePPM = scores$usePPM, useAdduct = scores$useAdduct,
                        brackets_ignore = scores$brackets_ignore)

    object <- labelRows(object, useID = labels$useID,minScore = labels$minScore,
                        maxRankX = labels$maxRankX, maxRankY = labels$maxRankY,
                        method = labels$method, balanced = labels$balanced,
                        remove = labels$remove, delta = labels$delta,
                        maxRTerr = labels$maxRTerr, rtOrder = labels$rtOrder,
                        resolveConflicts = labels$resolveConflicts,
                        brackets_ignore = labels$brackets_ignore)

    if(isTRUE(union))
        object <- updateTables(object, xdata = xdata, ydata = ydata)

    return(object)
}
