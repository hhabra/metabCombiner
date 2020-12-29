#' Format List of RT Ordered Pairs
#'
#' @param object metabCombiner object
#'
#' @param anchors data.frame of ordered pair features for RT fitting
#'
#' @param weights numeric fit weights
#'
#' @param useID  logical. Option to use matched IDs to inform fit
#'
#' @noRd
formatAnchors <- function(object, anchors, weights, useID){
    if(!useID | is.null(anchors[["labels"]]))
        anchors[["labels"]] = "A"
    cTable <- combinedTable(object)
    rtx <- c(min(cTable[["rtx"]]), anchors[["rtx"]], max(cTable[["rtx"]]))
    rty <- c(min(cTable[["rty"]]), anchors[["rty"]], max(cTable[["rty"]]))
    labels <- c("I", anchors[["labels"]], "I")
    if(length(weights) == nrow(anchors))
        weights <- c(1, weights, 1)

    rts <- data.frame(rtx, rty, labels, weights, stringsAsFactors = FALSE)
    return(rts)
}


#' Detect Outliers from Residuals
#'
#' @param residuals Matrix of numeric residuals
#'
#' @param include integer indices of non-outlier anchors
#'
#' @param vals numeric vector: k values for GAM fits, spans for loess fits
#'
#' @param outlier Thresholding method for outlier dection. If "MAD", the
#' threshold is the mean absolute deviation (MAD) times \code{coef}; if
#' "boxplot", the threshold is \code{coef} times IQR plus 3rd quartile of
#' a model's absolute residual values.
#'
#' @param coef numeric (> 1) multiplier for determining thresholds for outliers
#' (see \code{outlier} argument)
#'
#' @noRd
flagOutliers <- function(residuals, include, vals, outlier, coef)
{
    if(outlier == "MAD")
        thresholds <- coef * matrixStats::colMeans2(residuals, rows = include)
    else if(outlier == "boxplot")
        thresholds <- matrixStats::colIQRs(residuals, rows = include) * coef +
            matrixStats::colQuantiles(residuals, probs = 0.75, rows = include)
    flags <- vapply(seq(1,length(vals)), function(i){
        fl <- residuals[,i] >  thresholds[i]
        return(fl)
    }, logical(nrow(residuals)))

    return(flags)
}


#' @title Filter Outlier Ordered Pairs
#'
#' @description
#' Helper function for \code{\link{fit_gam}} & \code{\link{fit_loess}}. It
#' filters the set of ordered pairs using the residuals calculated from
#' multiple GAM / loess fits.
#'
#' @param rts  Data frame of ordered retention time pairs.
#'
#' @param fit  Either "gam" for GAM fits, or "loess" for loess fits
#'
#' @param vals numeric values: k values for GAM fits, spans for loess fits
#'
#' @param outlier Thresholding method for outlier dection. If "MAD", the
#' threshold is the mean absolute deviation (MAD) times \code{coef}; if
#' "boxplot", the threshold is \code{coef} times IQR plus 3rd quartile of
#' a model's absolute residual values.
#'
#' @param coef numeric (> 1) multiplier for determining thresholds for outliers
#' (see \code{outlier} argument)
#'
#' @param iterFilter integer number of outlier filtering iterations
#'
#' @param prop  numeric. A point is excluded if deemed a residual in more than
#' this proportion of fits. Must be between 0 & 1.
#'
#' @param bs character. Choice of spline method from mgcv; either "bs" or "ps"
#'
#' @param m integer. Basis and penalty order for GAM; see ?mgcv::s
#'
#' @param family character. Choice of mgcv family; see: ?mgcv::family.mgcv
#'
#' @param method character. Smoothing parameter estimation method; see:
#' ?mgcv::gam
#'
#' @param optimizer character. Method to optimize smoothing parameter; see:
#' ?mgcv::gam
#'
#' @param loess.pars parameters for LOESS fitting; see ?loess.control
#'
#' @param message Option to print message indicating function progress
#'
#' @param ... other arguments passed to \code{mgcv::gam}.
#'
#' @return anchor rts data frame with updated weights.
filterAnchors <- function(rts, fit, vals, outlier, coef, iterFilter,
                          prop, bs, m, family, method, optimizer, loess.pars,
                          message,...)
{
    iter <- 0
    while(iter < iterFilter){
        if(message)
            cat("Performing filtering iteration: ", iter + 1, "\n", sep = "")
        if(fit == "gam")
            residuals <- suppressWarnings(vapply(vals, function(v){
                model <- mgcv::gam(rty ~ s(rtx, k = v, bs = bs, m = m,...),
                                data = rts, family = family, method = method,
                                weights = rts[["weights"]],
                                optimizer = c("outer", optimizer), ...)

                res <- abs(model[["fitted.values"]] - model[["y"]])
                return(res)
            }, numeric(nrow(rts))))
        else if(fit == "loess")
            residuals <- suppressWarnings(vapply(vals, function(v){
                model <- stats::loess(rty ~ rtx, data = rts, span = v,
                                    degree = 1, weights = rts$weights,
                                    control = loess.pars, family = "s")

                res <- abs(model[["residuals"]])
                return(res)
            }, numeric(nrow(rts))))

        include <- which(rts[["weights"]] != 0)
        flags <- flagOutliers(residuals, include, vals, outlier, coef)
        remove <- rowSums(flags)/ncol(flags) > prop
        remove[rts[["labels"]] == "I"] <- FALSE

        if(sum(!remove) > 22 & sum(!remove) < length(remove))
            rts$weights[remove] <- 0
        else
            break
        iter <- iter+1
    }
    return(rts)
}

#' @title Cross Validation for Model Fits
#'
#' @description
#' Helper function for \code{fit_gam()} & \code{fit_loess()}. Determines
#' optimal value of \code{k} basis functions for Generalized Additive Model
#' fits or \code{span} for loess fits from among user-defined choices,
#' using a 10-fold cross validation minimizing mean squared error.
#'
#' @param rts data.frame of ordered pair retention times
#'
#' @param fit  Either "gam" for GAM fits, or "loess" for loess fits
#'
#' @param vals numeric vector: k values for GAM fits, spans for loess fits.
#' Best value chosen by 10-fold cross validation.
#'
#' @param bs character. Choice of spline method, either "bs" or "ps"
#'
#' @param family character. Choice of mgcv family; see: ?mgcv::family.mgcv
#'
#' @param m integer. Basis and penalty order for GAM; see ?mgcv::s
#'
#' @param method character. Smoothing parameter estimation method; see:
#' ?mgcv::gam
#'
#' @param optimizer  character. Method to optimize smoothing parameter; see:
#' ?mgcv::gam
#'
#' @param loess.pars parameters for LOESS fitting; see ?loess.control
#'
#' @param message Option to print message indicating function progress
#'
#' @param ... Other arguments passed to \code{mgcv::gam}.
#'
#' @return Optimal parameter value as determined by 10-fold cross validation
crossValFit <- function(rts, fit, vals, bs, family, m, method, optimizer,
                        loess.pars, message, ...)
{
    if(message) cat("Performing 10-fold cross validation\n")
    rts <- dplyr::filter(rts, .data$weights != 0)
    N <- nrow(rts) - 1
    folds <- caret::createFolds(seq(2,N), k = 10, returnTrain = FALSE)

    cv_errors <- vapply(folds, function(f){
        f <- f + 1
        rts_train <- rts[-f,]
        rts_test <- rts[f,]

        #error for each span for fold f
        errors <- suppressWarnings(vapply(vals, function(v){
            if(fit == "gam")
                model <- mgcv::gam(rty ~ s(rtx, k = v, bs = bs, m = m, ...),
                            data = rts_train, family = family, method = method,
                            optimizer = c("outer",optimizer),
                            weights = rts_train$weights, ...)
            else if (fit == "loess")
                model <- stats::loess(rty ~ rtx, data = rts_train, span = v,
                                    degree = 1, weights = rts_train$weights,
                                    control = loess.pars,  family = "s")

            preds <- stats::predict(model, newdata = rts_test)
            MSE <- sum((preds - rts_test[["rty"]])^2)/ length(preds)

            return(MSE)
        }, numeric(1)))

        return(errors)
    }, numeric(length(vals)))

    mean_cv_errors <- matrixStats::rowMeans2(as.matrix(cv_errors))
    best_val <- vals[which.min(mean_cv_errors)]
    return(best_val)
}

#' @title Fit RT Projection Model With GAMs
#'
#' @description
#' Fits a (penalized) basis splines curve through a set of ordered pair
#' retention times, modeling one set of retention times (rty) as a function
#' on the other set (rtx). Outlier filtering iterations are performed first,
#' then with the remaining points, the best value of parameter \code{k} is
#' selected through 10-fold cross validation.
#'
#' @param object  a \code{metabCombiner} object.
#'
#' @param useID  logical. If set to TRUE, matched ID anchors detected from
#' previous step will never be flagged as outliers.
#'
#' @param k  integer k values controlling the dimension of the basis of the
#' GAM fit (see: ?mgcv::s). Best value chosen by 10-fold cross validation.
#'
#' @param iterFilter integer number of outlier filtering iterations to perform
#'
#' @param outlier Thresholding method for outlier dection. If "MAD", the
#' threshold is the mean absolute deviation (MAD) times \code{coef}; if
#' "boxplot", the threshold is \code{coef} times IQR plus 3rd quartile of
#' a model's absolute residual values.
#'
#' @param coef numeric (> 1) multiplier for determining thresholds for outliers
#' (see \code{outlier} argument)
#'
#' @param prop  numeric. A point is excluded if deemed a residual in more than
#' this proportion of fits. Must be between 0 & 1.
#'
#' @param bs   character. Choice of spline method from mgcv, either "bs" (basis
#' splines) or "ps" (penalized basis splines)
#'
#' @param family character. Choice of mgcv family; see: ?mgcv::family.mgcv
#'
#' @param weights Optional user supplied weights for each ordered pair. Must be
#' of length equal to number of anchors (n) or a divisor of (n + 2).
#'
#' @param m  integer. Basis and penalty order for GAM; see ?mgcv::s
#'
#' @param method character smoothing parameter estimation method; see:
#' ?mgcv::gam
#'
#' @param optimizer character. Method to optimize smoothing parameter; see:
#' ?mgcv::gam
#'
#' @param message Option to print message indicating function progress
#'
#' @param ... Other arguments passed to \code{mgcv::gam}.
#'
#' @details
#' A set of ordered pair retention times must be previously computed using
#' \code{selectAnchors()}. The minimum and maximum retention times from both
#' input datasets are included in the set as ordered pairs (min_rtx, min_rty)
#' & (max_rtx, max_rty). The \code{weights} argument initially determines the
#' contribution of each point to the model fits; they are equally weighed by
#' default, but can be changed using an \code{n+2} length vector, where n is
#' the number of ordered pairs and the first and last of the weights determines
#' the contribution of the min and max ordered pairs; by default, all weights
#' are initially set to 1 for equal contribution of each point.
#'
#' The model complexity is determined by \code{k}. Multiple values of k are
#' allowed, with the best value chosen by 10 fold cross validation. Before
#' this happens, certain ordered pairs are removed based on the model errors.
#' In each iteration, a GAM is fit using each selected value of k. Depending on
#' the \code{outlier} argument, a point is "removed" from the model (i.e. its
#' corresponding weight set to 0) if its residual is above the threshold
#' for a proportion of fitted models, as determined by \code{prop}. If an anchor
#' is an "identity" (idx = idy, detected in the \code{selectAnchors} by setting
#' \code{useID} to TRUE), then setting \code{useID} here prevents its removal.
#'
#' Other arguments, e.g. \code{family}, \code{m}, \code{optimizer}, \code{bs},
#' and \code{method} are GAM specific parameters from the \code{mgcv} R package.
#' The \code{family} option is currently limited to the "scat" (scaled t) and
#' "gaussian" families; scat family model fits are more robust to outliers than
#' gaussian fits, but compute much slower. Type of splines are currently limited
#' to basis splines ("bs" or "ps").
#'
#' @return metabCombiner with a fitted GAM model object
#'
#' @seealso
#' \code{\link{selectAnchors}},\code{\link{fit_loess}},
#'
#' @examples
#' data(plasma30)
#' data(plasma20)
#'
#' p30 <- metabData(plasma30, samples = "CHEAR")
#' p20 <- metabData(plasma20, samples = "Red", rtmax = 17.25)
#' p.comb = metabCombiner(xdata = p30, ydata = p20, binGap = 0.0075)
#'
#' p.comb = selectAnchors(p.comb, tolmz = 0.003, tolQ = 0.3, windy = 0.02)
#' anchors = getAnchors(p.comb)
#'
#' #version 1: using faster, but less robust, gaussian family
#' p.comb = fit_gam(p.comb, k = c(10,12,15,17,20), prop = 0.5,
#'     family = "gaussian", outlier = "MAD", coef = 2)
#'
#' \donttest{
#' #version 2: using slower, but more robust, scat family
#' p.comb = fit_gam(p.comb, k = seq(12,20,2), family = "scat",
#'                      iterFilter = 1, coef = 3, method = "GCV.Cp")
#'
#' #version 3 (with identities)
#' p.comb = selectAnchors(p.comb, useID = TRUE)
#' anchors = getAnchors(p.comb)
#' p.comb = fit_gam(p.comb, useID = TRUE, k = seq(12,20,2), iterFilter = 1)
#'
#' #version 4 (using identities and weights)
#' weights = ifelse(anchors$labels == "I", 2, 1)
#' p.comb = fit_gam(p.comb, useID = TRUE, k = seq(12,20,2),
#'                      iterFilter = 1, weights = weights)
#'
#' #version 5 (using boxplot-based outlier detection
#' p.comb = fit_gam(p.comb, k = seq(12,20,2), outlier = "boxplot", coef = 1.5)
#'
#' #to preview result of fit_gam
#' plot(p.comb, pch = 19, outlier = "h", xlab = "CHEAR Plasma (30 min)",
#'      ylab = "Red-Cross Plasma (20 min)", main = "Example GAM Fit")
#' }
#'
#' @export
fit_gam <- function(object, useID = FALSE, k = seq(10,20, by = 2),
                    iterFilter = 2, outlier = c("MAD", "boxplot"), coef = 2,
                    prop = 0.5, weights = 1, bs = c("bs", "ps"),
                    family = c("scat", "gaussian"), m = c(3,2),
                    method = "REML", optimizer = "newton", message = TRUE, ...)
{
    combinerCheck(isMetabCombiner(object), "metabCombiner")
    anchors <- getAnchors(object)
    check_fit_pars(anchors = anchors, fit = "gam", useID = useID, k = k,
                    iterFilter = iterFilter, coef = coef, prop = prop)
    outlier <- match.arg(outlier)
    bs <- match.arg(bs)
    family <- match.arg(family)

    rts <- formatAnchors(object, anchors, weights, useID)
    rts <- filterAnchors(rts = rts, fit = "gam", vals = k, outlier = outlier,
                         coef = coef, iterFilter = iterFilter, prop = prop,
                         bs = bs, m = m, family = family, method = method,
                         optimizer = optimizer, message = message, ...)

    if(length(k) > 1)
        best_k <- crossValFit(rts = rts, vals = k, fit = "gam", bs = bs, m = m,
                                family = family, method = method,
                                optimizer = optimizer, message = message,...)
    else
        best_k <- k

    if(message) cat("Fitting Model with k =", best_k, "\n")
    best_model <- mgcv::gam(rty ~ s(rtx, k = best_k, bs = bs, m = m, ...),
                            data = rts, family = family, method = method,
                            optimizer = c("outer", optimizer),
                            weights = rts[["weights"]], ...)

    anchors[["rtProj"]] <- stats::predict(best_model, newdata = anchors)
    object <- update_mc(object, anchors = anchors, fit = "gam",
                        model = best_model,stats = "best_k", values = best_k)
    return(object)
}


#' @title Fit RT Projection Model With LOESS
#'
#' @description
#' Fits a local regression smoothing spline through a set of ordered pair
#' retention times. modeling one set of retention times (rty) as a function
#' on the other set (rtx). Filtering iterations of high residual points are
#' performed first. Multiple acceptable values of \code{span} can be used, with
#' one value selected through 10-fold cross validation.
#'
#' @param object  a \code{metabCombiner} object.
#'
#' @param useID  logical. If set to TRUE, matched ID anchors detected from
#' previous step will never be flagged outliers.
#'
#' @param spans numeric span values (between 0 & 1) used for loess fits
#'
##' @param outlier Thresholding method for outlier dection. If "MAD", the
#' threshold is the mean absolute deviation (MAD) times \code{coef}; if
#' "boxplot", the threshold is \code{coef} times IQR plus 3rd quartile of
#' a model's absolute residual values.
#'
#' @param coef numeric (> 1) multiplier for determining thresholds for outliers
#' (see \code{outlier} argument)
#'
#' @param iterFilter integer number of outlier filtering iterations to perform
#'
#' @param prop  numeric. A point is excluded if deemed a residual in more than
#' this proportion of fits. Must be between 0 & 1.
#'
#' @param iterLoess  integer. Number of robustness iterations to perform in
#'                   \code{loess()}.See ?loess.control for more details.
#'
#' @param weights Optional user supplied weights for each ordered pair. Must be
#' of length equal to number of anchors (n) or a divisor of (n + 2)
#'
#' @param message Option to print message indicating function progress
#'
#' @return \code{metabCombiner} object with \code{model} slot updated to
#' contain the fitted loess model
#'
#' @seealso
#' \code{\link{selectAnchors}},\code{\link{fit_gam}}
#'
#' @examples
#' data(plasma30)
#' data(plasma20)
#'
#' p30 <- metabData(plasma30, samples = "CHEAR")
#' p20 <- metabData(plasma20, samples = "Red", rtmax = 17.25)
#' p.comb = metabCombiner(xdata = p30, ydata = p20, binGap = 0.0075)
#' p.comb = selectAnchors(p.comb, tolmz = 0.003, tolQ = 0.3, windy = 0.02)
#'
#' #version 1
#' p.comb = fit_loess(p.comb, spans = seq(0.2,0.3,0.02), iterFilter = 1)
#'
#' #version 2 (using weights)
#' anchors = getAnchors(p.comb)
#' weights = c(2, rep(1, nrow(anchors)), 2)  #weight = 2 to boundary points
#' p.comb = fit_loess(p.comb, spans = seq(0.2,0.3,0.02), weights = weights)
#'
#' #version 3 (using identities)
#' p.comb = selectAnchors(p.comb, useID = TRUE, tolmz = 0.003)
#' p.comb = fit_loess(p.comb, spans = seq(0.2,0.3,0.02), useID = TRUE)
#'
#' #to preview result of fit_loess
#' plot(p.comb, fit = "loess", xlab = "CHEAR Plasma (30 min)",
#'      ylab = "Red-Cross Plasma (20 min)", pch = 19,
#'      main = "Example fit_loess Result Fit")
#'
#' @export
fit_loess <- function(object, useID = FALSE, spans = seq(0.2, 0.3, by = 0.02),
                    outlier = c("MAD", "boxplot"), coef = 2, iterFilter = 2,
                    prop = 0.5, iterLoess = 10, weights = 1, message = TRUE)
{
    combinerCheck(isMetabCombiner(object), "metabCombiner")
    anchors <- object@anchors

    check_fit_pars(anchors = anchors, fit = "loess", useID = useID,
                    iterFilter = iterFilter, coef = coef, prop = prop,
                    iterLoess = iterLoess, spans = spans)

    outlier <- match.arg(outlier)
    loess.pars <- loess.control(iters = iterLoess, surface = "direct")
    rts <- formatAnchors(object, anchors, weights, useID)

    rts <- filterAnchors(rts = rts, fit = "loess", iterFilter = iterFilter,
                        outlier = outlier, coef = coef, prop = prop,
                        vals = spans, loess.pars = loess.pars,
                        message = message)

    if(length(spans) > 1)
        best_span <- crossValFit(rts = rts, fit = "loess", vals = spans,
                                loess.pars = loess.pars, message = message)
    else
        best_span <- spans

    cat("Fitting Model with span =", best_span,"\n")

    best_model <- loess(rty ~ rtx, data = rts, span = best_span, degree = 1,
                        family = "symmetric", control = loess.pars,
                        weights = rts[["weights"]])

    anchors[["rtProj"]] <- stats::predict(best_model, anchors)
    object <- update_mc(object, anchors = anchors, fit = "loess",
                        model = best_model, stats = "best_span",
                        values = best_span)
    return(object)
}




