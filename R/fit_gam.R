##
#' @title Filter Outlier Ordered Pairs
#'
#' @description
#' This is a helper function for \code{fit_gam}. It filters retention time ordered
#' pairs using the residuals calculated from multiple GAM fits.
#'
#' @param rts  Data frame of ordered retention time pairs.
#'
#' @param k  Integer vector. Selection of values to control the number of basis
#' functions for GAM construction.
#'
#' @param bs   character. Choice of spline method from mgcv; either "bs" or "ps".
#'
#' @param iterFilter integer. Number of residual filtering iterations to perform.
#'
#' @param family character. Choice of mgcv family; see: ?mgcv::family.mgcv
#'
#' @param m integer. Basis and penalty order for GAM; see ?mgcv::s
#'
#' @param method character. Smoothing parameter estimation method; see:
#' ?mgcv::gam
#'
#' @param optimizer character. Method to optimize smoothing parameter; see:
#' ?mgcv::gam
#'
#' @param ratio numeric. Residual multiplier for determining outliers.
#'
#' @param frac  numeric. A point is excluded if deemed a residual in more than
#' this fraction value times the number of fits. Must be between 0 & 1.
#'
#' @param ... Other arguments passed to \code{mgcv::gam}.
#'
#' @return
#' Ordered retention time pairs data frame with updated weights.
#'
##
filterAnchorsGAM <- function(rts, k, bs, iterFilter, family, m, method, optimizer,
                             ratio, frac,...){
    iteration = 0

    while(iteration < iterFilter){
        cat("Performing filtering iteration: ", iteration + 1, "\n", sep = "")
        residuals <- sapply(k, function(ki){
              model <- mgcv::gam(rty ~ s(rtx, k = ki, bs = bs, m = m,...),
                             data = rts, family = family, method = method,
                             weights = rts[["weights"]],
                             optimizer = c("outer", optimizer), ...)

              res = abs(model[["fitted.values"]] - model[["y"]])
              return(res)
        })

        include = which(rts[["weights"]] == 1)

        ##flagging excessively high residuals
        thresholds <- ratio * colMeans(residuals[include,] %>% as.matrix())

        flags <- sapply(1:length(k), function(i){
            fl <- (residuals[,i] > thresholds[i])
            return(fl)
        })

        ##remove features with high proportion of excessive residuals
        fracs = rowSums(flags)/ncol(flags)
        remove = fracs > frac
        remove[rts[["labels"]] == "I"] = FALSE

        ##assigning 0 weight to outlier anchors;
        ##must retain 20 anchors plus both endpoints
        if(sum(!remove) > 22 & sum(!remove) < length(remove))
            rts$weights[remove] = 0
        else
            break

        iteration = iteration+1
    }

    return(rts)
}

##
#' @title Cross Validation for GAM model
#'
#' @description
#' Helper function for \code{fit_gam()}. Determines optimal value
#' of \code{k} basis functions (among user-defined choices) for Generalized
#' Additive Model, using a 10-fold cross validation.
#'
#' @param rts Data frame of ordered retention time pairs.
#'
#' @param k  Integer vector. Selection of values to control the number of basis
#' functions for gam construction.
#'
#' @param bs   character. Choice of spline method from mgcv, either "bs" or "ps"
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
#' @param seed integer. Psuedo-random seed generator for cross-validation.
#'
#' @param ... Other arguments passed to \code{mgcv::gam}.
#'
#' @return Optimal value for \code{k} as determined by 10-fold cross validation.
#'
##
crossValidationGAM <- function(rts, k, bs, m, family, method, optimizer, seed,...)
{
    cat("Performing 10-fold cross validation.")

    #excluding filtered points
    rts = dplyr::filter(rts, weights != 0)

    ##skipping last entry
    N = nrow(rts) - 1

    #reproducibility seed
    set.seed(seed)

    #creating 10 folds; boundary ordered pairs are skipped
    folds <- caret::createFolds(2:N, k = 10, returnTrain = FALSE)

    ##collecting cross validation errors
    cv_errors <- sapply(folds, function(f){
        f = f + 1
        rts_train <- rts[-f,]
        rts_test <- rts[f,]

        #error for each span for fold f
        errors <- sapply(k, function(ki){
            model <- mgcv::gam(rty ~ s(rtx, k = ki, bs = bs, m = m, ...),
                           data = rts_train, family = family, method = method,
                           optimizer = c("outer",optimizer),
                           weights = rts_train$weights,
                           ...)

            preds <- stats::predict(model, newdata = rts_test)

            MSE = sum((preds - rts_test[["rty"]])^2)/ length(preds)

            return(MSE)
        })

        return(errors)
    })

    mean_cv_errors = rowMeans(matrix(cv_errors, nrow = length(k)))

    best_k = k[which.min(mean_cv_errors)]

    return(best_k)
}

##
#' @title Fit RT Projection Model With GAMs
#'
#' @description
#' Fits a (penalized) basis splines curve through a set of ordered pair
#' retention times. One dataset's retention times (ydata) serve as the dependent
#' variable, fitted as a function of retention times of the other (xdata).
#'
#' Filtering iterations of high residual points are performed, controlled by
#' \code{iterFilter}, and multiple acceptable values of \code{k} used, with
#' one value selected using 10-fold cross validation.
#'
#' @param object  a metabCombiner object.
#'
#' @param useID  logical. Option to use matched IDs to inform modeling process.
#'
#' @param k   integer vector. Selection of values to control the number of basis
#' functions for gam construction. Best value chosen by 10-fold cross validation.
#'
#' @param ratio numeric. A point is considered an outlier if ratio of residual to
#' mean residual of a fit exceeds this value. Must be greater than 1.
#'
#' @param frac  numeric. A point is excluded if deemed a residual in more than
#' this fraction value times the number of fits. Must be between 0 & 1.
#'
#' @param iterFilter integer. Number of residual filtering iterations to perform.
#'
#' @param bs   character. Choice of spline method from mgcv, either "bs" (basis
#' splines) or "ps" (penalized basis splines)
#'
#' @param family character. Choice of mgcv family; see: ?mgcv::family.mgcv
#'
#' @param weights Numeric prior weights determines the contribution of each point
#' to the model. see: ?mgcv::gam
#'
#' @param m  integer. Basis and penalty order for GAM; see ?mgcv::s
#'
#' @param method  character. Smoothing parameter estimation method; see:
#' ?mgcv::gam
#'
#' @param optimizer character. Method to optimize smoothing parameter; see:
#' ?mgcv::gam
#'
#' @param seed integer. Psuedo-random seed generator for cross validation.
#'
#' @param ... Other arguments passed to \code{mgcv::gam}.
#'
#' @details
#' A set of ordered pair retention times must be previously computed using
#' \code{selectAnchors()}. The minimum and maximum retention times from both
#' input datasets are included in the set (min_rtx, min_rty) & (max_rtx, max_rty).
#'
#' The \code{weights} argument initially determines the contribution of each point
#' to the model fits; they are equally weighed by default, but can be changed using
#' an \code{n+2} length vector, where n is the number of ordered pairs and the
#' first and last of the weights determines the contribution of the min and max
#' ordered pairs.
#'
#' The model complexity is determined by argument \code{k}. Multiple values of k
#' are expected, with the best value chosen by 10 fold cross validation. Before
#' this happens, certain ordered pairs are removed based on the model errors.
#' In each iteration, a GAM is fit using each selected value of k. A point is
#' "removed" (its corresponding \code{weights} value set to 0) if its residual
#' is \code{ratio} times average residual for a fraction of fitted models, as
#' determined by \code{frac}. If an ordered pair is an "identity" (discovered
#' in the \code{selectAnchors} by setting the \code{useID} to TRUE), then
#' setting \code{useID} here will prevent its removal.
#'
#' Other arguments, such as \code{family}, \code{m}, \code{optimizer}, \code{bs},
#' and \code{method} are GAM specific parameters. Type of splines are currently
#' limited to basis splines (\eqn{bs = "bs"}) or penalized basis splines
#' (\eqn{bs = "ps"}).
#'
#' @return metabCombiner object with a fitted GAM model
#'
#' @seealso
#' [crossValidationGAM], [filterAnchorsGAM], [fit_loess]
#'
#' @export
##
fit_gam <- function(object, useID = FALSE, k = seq(10,20, by = 2), iterFilter = 2,
                    ratio = 2, frac = 0.5, bs = c("bs", "ps"),
                    family = c("scat", "gaussian"), weights = 1, m = c(3,2),
                    method = "REML", optimizer = "newton", seed = 100, ...)
{
    code = isMetabCombiner(object)

    if(code)
        stop(combinerError(code, "metabCombiner"))

    anchors = object@anchors

    ##parameter value checks
    if(nrow(anchors) == 0)
        stop("Missing 'anchors' in metabCombiner object. selectAnchors() has not
              been called yet")

    if(nrow(anchors) <= 20)
        stop("Number of anchors too low (less than 20). Examine parameters from
              selectAnchors().")

    if(class(useID) != "logical")
        stop("parameter 'useID' must be logical")

    if(!is.numeric(k) | any(k < 3 | k >= nrow(anchors)))
        stop("all values of k must be between 3 and number of anchors")

    if(class(iterFilter) != "numeric" | as.integer(iterFilter) < 0){
        warning("parameter 'iterFilter' must be positive. Setting value to
               default(2)")
        iterFilter = 2
    }

    if(class(ratio)!= "numeric" | ratio < 1){
        warning("Invalid parameter 'ratio'. Must be a number greater than 1.
                Setting value to default(2).")
        ratio = 2
    }

    if(frac <= 0 | frac > 1 | !is.numeric(frac)){
        warning("parameter 'frac' must be numeric between 0 & 1. Setting value
                to default(0.5)")
        frac = 0.5
    }

    if(useID != TRUE)
        anchors[["labels"]] = "A"

    ##gam parameters
    bs = match.arg(bs)

    family = match.arg(family)
    if(family == "scat") family = scat()
    else if (family == "gaussian") family = "gaussian"

    ##appending minimum and maximum RT values to retention time lists
    cTable = combinerTable(object)
    RTx = c(base::min(cTable[["rtx"]]), base::max(cTable[["rtx"]]))
    RTy = c(base::min(cTable[["rty"]]), base::max(cTable[["rty"]]))
    rtx <- c(RTx[1], anchors[["rtx"]], RTx[2])
    rty <- c(RTy[1], anchors[["rty"]], RTy[2])
    labels <- c("I", anchors[["labels"]], "I")

    rts <- data.frame(rtx, rty, labels, weights = weights,
                      stringsAsFactors = FALSE)

    ###Anchor filtering procedure
    rts = suppressWarnings(filterAnchorsGAM(rts = rts,
                                            k = k,
                                            bs = bs,
                                            m = m,
                                            family = family,
                                            method = method,
                                            optimizer = optimizer,
                                            iterFilter = iterFilter,
                                            ratio = ratio,
                                            frac = frac,...))

    ##10 fold cross validation: find the best value for k
    best_k <- suppressWarnings(crossValidationGAM(rts = rts,
                                                  k = k,
                                                  bs = bs,
                                                  m = m,
                                                  family = family,
                                                  method = method,
                                                  optimizer = optimizer,
                                                  seed = seed,...))



    best_model <- mgcv::gam(rty ~ s(rtx, k = best_k, bs = bs, m = m, ...),
                                  data = rts, family = family, method = method,
                                  optimizer = c("outer", optimizer),
                                  weights = rts[["weights"]], ...)

    object@model[["gam"]] = best_model
    object@stats[["best_k"]] = best_k

    return(object)
}

