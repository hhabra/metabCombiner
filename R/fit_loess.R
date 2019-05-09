##filterAnchors
#'
#' This is a helper function for the fit_loess function.
#' Filters anchor retention time ordered pairs using the residuals calculated
#' from an initial set of loess fits (using specified spans). Those residuals
#' that exceed some ratio * average_residuals some
#'
##
filterAnchors <- function(rts, spans, iterFilter, ratio, frac, loess.parameters){
    iteration = 0
    
    while(iteration < iterFilter){
        residuals <- lapply(spans, function(s){
              model <- loess(rty ~ rtx, data = rts, span = s, degree = 1,
                     family = "symmetric", control = loess.parameters)
      
              res = model$residuals
              return(res)
        })
    
        residuals = matrix(unlist(residuals), ncol = length(spans)) %>% abs()
    
        ##flagging excessively high residuals
        thresholds <- colMeans(residuals) * ratio
    
        flags <- lapply(1:length(spans), function(i){
            fl <- (residuals[,i] > thresholds[i])
            return(fl)
        }) %>% unlist %>% matrix(ncol = length(spans))
    
        ##remove features with high proportion of excessive residuals
        fracs = rowSums(flags)/ncol(flags)
        keep = fracs < frac
        keep[rts$labels == "I"] = TRUE
    
        ##filtering ordered pairs; need 20 anchors plus both endpoints
        if(sum(keep) > 22 & sum(keep) < length(keep))          
            rts = rts[keep,]
         else 
            break
        
        iteration = iteration+1
    }
  
  return(rts)
}

##crossValidation
#'
#'This is a helper function for the fit_loess function.
#'Performs cross validation to choose the most appropriate span loess for loess fit. 
#'
#'
##

crossValidationLoess <- function(rts, spans, seed, loess.parameters){
    
    ##skipping last entry 
    N = nrow(rts) - 1
  
    #reproducibility seed                                                        
    set.seed(seed)
  
    #creating 10 folds; boundary ordered pairs are skipped
    folds <- caret::createFolds(2:N, k = 10, returnTrain = FALSE)
  
    ##collecting cross validation errors
    cv_errors <- lapply(folds, function(f){
        f = f + 1      
        rts_train <- rts[-f,]
        rts_test <- rts[f,]
    
        #error for each span for fold f
        errors <- sapply(spans, function(s){
            model <- loess(rty ~ rtx, data = rts_train, span = s, degree = 1,
                     family = "symmetric",control = loess.parameters)
      
            preds <- predict(model, newdata = rts_test)
      
            MSE = sum((preds - rts_test$rty)^2)/ length(preds)
      
            return(MSE)
        })
    
        return(errors)
    }) %>% as.data.frame()
  
    return(cv_errors)
}

##fit_loess
#'
#' @param object  metabCombiner object with non-empty anchors field.
#' 
#' @param useID   logical. Option to not filter identity-labeled ordered pair 
#'
#' @param spans   numeric vector. Span values (between 0 & 1) used for loess fits.
#'
#' @param ratio   numeric. Residual multiplier for determining outliers.
#' 
#' @param frac    numeric. Threshold proportion of fits needed to determine if 
#'                ordered pair is an outlier.
#'                
#' @param iterFilter  integer. Number of residual filtering iterations to perform. 
#' 
#' @param iterLoess  integer. Number of robustness iterations to perform in loess().
#'                   See ?loess.control for more details.
#'                   
#' @param seed    integer. Psuedo-random seed generator for cross validation folds. 
#'
##
fit_loess <- function(object, useID = FALSE, spans = seq(0.2, 0.4, by = 0.05), 
                      ratio = 2, frac = 0.5, iterFilter = 2, iterLoess = 10, 
                      seed = 100)
{
    if(class(object) != "metabCombiner")
        stop(base::paste(object, "is not a metabCombiner object"), sep = " ")
  
    anchors = object@anchors
    
    ##parameter value checks 
    if(nrow(anchors) == 0)
        stop("Missing 'anchors' in metabCombiner object. selectAnchors() has not 
              been called yet")
  
    if(nrow(anchors) <= 20)
        stop("Number of anchors too low (less than 20). Examine parameters from 
              selectAnchors().")
  
    if(any(spans <= 0 | spans >= 1) | !is.numeric(spans))
        stop("all values in parameter 'spans' must be numeric between 0 & 1")
    
    if(class(useID) != "logical")
        stop("parameter 'useID' must be logical")
    
    if(class(ratio)!= "numeric" | ratio < 1){
        warning("Invalid parameter 'ratio'. Must be a number greater than 1. 
                Setting value to default(2).")
        ratio = 2
    }
  
    if(class(iterFilter) != "numeric" | as.integer(iterFilter) < 0){
        warning("parameter 'iterFilter' must be positive. Setting value to 
              default(2)")
        iterFilter = 2
    }
    
    if(class(iterLoess) != "numeric" | as.integer(iterLoess) < 1){
        warning("parameter 'iterLoess' must be at least 1. Setting value to 
                default(10)")
        iterLoess = 10
    }
  
    if(frac <= 0 | frac > 1 | !is.numeric(frac)){
        warning("parameter 'frac' must be numeric between 0 & 1. Setting value to
                 default(0.5)")
        frac = 0.5
    }
    
    if(!useID)
        anchors$labels = "A"
  
    loess.parameters = stats::loess.control(iterations = iterLoess, surface = "direct")
  
    ##appending minimum and maximum RT values to retention time lists
    
    cTable = combinerTable(object)
    minRTx = base::min(cTable$rtx)
    minRTy = base::min(cTable$rty)
    maxRTx = base::max(cTable$rtx)
    maxRTy = base::max(cTable$rty)
    rtx <- c(minRTx, anchors$rtx, maxRTx)
    rty <- c(minRTy, anchors$rty, maxRTy)
    labels <- c("I", anchors$labels, "I")
    rts <- data.frame(rtx, rty, labels, stringsAsFactors = FALSE)
  
    ###Anchor filtering procedure
    rts = suppressWarnings(filterAnchors(rts, spans, iterFilter = iterFilter, 
                                         ratio, frac, loess.parameters))
  
    ##10 fold cross validation: find the best span, s
    cv_errors <- suppressWarnings(crossValidationLoess(rts, spans, seed, loess.parameters))
    mean_cv_errors = rowMeans(cv_errors)
  
    best_span = spans[which.min(mean_cv_errors)]
  
    best_model <- loess(rty ~ rtx, data = rts, span = best_span, degree = 1,
                      family = "symmetric",control = loess.parameters)
  
    object@model$loess = best_model
  
    return(object)
}


