## Contains a LOESS fitting method

#' @title Filter Outlier Ordered Pairs
#'
#' @description
#' This is a helper function for \code{fit_loess}. It filters retention time
#' ordered pairs using the residuals calculated from multiple LOESS fits.
#'
#' @param rts Data frame of ordered retention time pairs.
#'
#' @param spans numeric vector. Selection of values (between 0 & 1) to control
#' the degree of smoothing of loess curve.
#'
#' @param iterFilter integer. Number of residual filtering iterations to perform.
#'
#' @param ratio   numeric. A point is considered an outlier if ratio of residual to
#' mean residual of a fit exceeds this value. Must be greater than 1.
#'
#' @param frac   numeric. A point is excluded if deemed a residual in more than
#' this fraction value times the number of fits. Must be between 0 & 1.
#'
#' @param loess.parameters  Parameters for LOESS fitting. See ?loess.control.
#'
#' @return
#' Ordered retention time pairs data frame with updated weights.
#'
#' @seealso
#' [filterAnchorsGAM],[fit_loess], [crossValidationLoess]
##
filterAnchorsLoess <- function(rts, spans, iterFilter, ratio, frac,loess.parameters)
{
    iteration = 0

    while(iteration < iterFilter){
        cat("Performing filtering iteration: ", iteration + 1, "\n", sep = "")

        residuals <- sapply(spans, function(s){
            model <- stats::loess(rty ~ rtx, data = rts, span = s, degree = 1,
                     family = "symmetric", control = loess.parameters,
                     weights = rts[["weights"]])

            res = abs(model[["residuals"]])
            return(res)
        })

        include = which(rts[["weights"]] == 1)

        ##flagging excessively high residuals
        thresholds <- ratio * colMeans(residuals[include,] %>% as.matrix())

        flags <- sapply(1:length(spans), function(i){
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

#' @title  Cross Validation for Loess Model
#'
#' @param rts Data frame of ordered retention time pairs.
#'
#' @param spans numeric vector. Selection of values (between 0 & 1) to control
#' the degree of smoothing of loess curve.
#'
#' @param seed integer. Psuedo-random seed generator for cross validation.
#'
#' @param loess.parameters Parameters for LOESS fitting. See ?loess.control.
#'
#' @return
#' Optimal value among pre-selected values for \code{spans}.
#'
#'
crossValidationLoess <- function(rts, spans, seed,loess.parameters){
    cat("Performing 10-fold cross validation.")

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
        errors <- sapply(spans, function(s){
            model <- stats::loess(rty ~ rtx, data = rts_train, span = s,
                           degree = 1, family = "symmetric",
                           control = loess.parameters,
                           weights = rts_train[["weights"]])

            preds <- stats::predict(model, newdata = rts_test)

            MSE = sum((preds - rts_test[["rty"]])^2)/ length(preds)

            return(MSE)
        })

        return(errors)
    })

    return(cv_errors)
}

#' @title Fit RT Projection Model With LOESS
#'
#' @description
#' Fits a local regression smoothing spline curve through a set of ordered pair
#' retention times. One dataset's retention times (ydata) serve as the dependent
#' variable, fitted as a function of retention times of the other (xdata).
#'
#' Filtering iterations of high residual points are performed, controlled by
#' \code{iterFilter}, and multiple acceptable values of \code{span} used, with
#' one value chosen using 10-fold cross validation.
#'
#' @param object  metabCombiner object.
#'
#' @param useID   logical. Option to not filter identity-labeled ordered pair.
#'
#' @param spans   numeric vector. Span values (between 0 & 1) used for loess fits
#'
#' @param ratio   numeric. A point is considered an outlier if ratio of residual to
#' mean residual of a fit exceeds this value. Must be greater than 1.
#'
#' @param frac   numeric. A point is excluded if deemed a residual in more than
#' this fraction value times the number of fits. Must be between 0 & 1.
#'
#' @param iterFilter  integer. Number of residual filtering iterations to perform.
#'
#' @param iterLoess  integer. Number of robustness iterations to perform in
#'                   \code{loess()}.See ?loess.control for more details.
#'
#' @param weights  numeric. Optional user supplied weights. Note: vector length
#'                  must be #num_anchors + 2 (first and last are endpoints).
#'
#' @param seed  integer. Psuedo-random seed generator for cross validation.
#'
#'
#' @export
fit_loess <- function(object, useID = FALSE, spans = seq(0.2, 0.4, by = 0.05),
                      iterFilter = 2, ratio = 2, frac = 0.5, iterLoess = 10,
                      weights = 1, seed = 100)
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

    if(any(spans <= 0 | spans >= 1) | !is.numeric(spans))
        stop("all values in parameter 'spans' must be numeric between 0 & 1")

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
       warning("parameter 'frac' must be numeric between 0 & 1. Setting value to
                 default(0.5)")
        frac = 0.5
    }

    if(class(iterLoess) != "numeric" | as.integer(iterLoess) < 1){
        warning("parameter 'iterLoess' must be at least 1. Setting value to
                default(10)")
        iterLoess = 10
    }

    if(!useID | is.null(anchors[["labels"]]))
        anchors[["labels"]] = "A"

    loess.parameters = loess.control(iterations = iterLoess, surface = "direct")

    ##appending minimum and maximum RT values to retention time lists
    cTable = combinerTable(object)

    RTx = c(base::min(cTable[["rtx"]]), base::max(cTable[["rtx"]]))
    RTy = c(base::min(cTable[["rty"]]), base::max(cTable[["rty"]]))
    rtx <- c(RTx[1], anchors[["rtx"]], RTx[2])
    rty <- c(RTy[1], anchors[["rty"]], RTy[2])
    labels <- c("I", anchors[["labels"]], "I")
    rts <- data.frame(rtx, rty, labels, weights = weights, stringsAsFactors = FALSE)

    #Anchor filtering procedure
    rts = suppressWarnings(filterAnchorsLoess(rts, spans, iterFilter = iterFilter,
                                         ratio, frac, loess.parameters))

    ##10 fold cross validation: find the best span, s
    cv_errors <- suppressWarnings(crossValidationLoess(rts, spans, seed,
                                                       loess.parameters))

    mean_cv_errors = rowMeans(matrix(cv_errors, nrow = length(spans)))
    best_span = spans[which.min(mean_cv_errors)]

    best_model <- loess(rty ~ rtx, data = rts, span = best_span, degree = 1,
                        family = "symmetric",control = loess.parameters,
                        weights = rts[["weights"]])

    object@model[["loess"]] = best_model
    object@stats[["best_span"]] = best_span

    return(object)
}


