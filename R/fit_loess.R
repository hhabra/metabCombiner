#'This is a helper function for the fit_loess function.
#' Filters anchor retention time ordered pairs using the residuals calculated
#' from an initial set of loess fits (using specified spans). Those residuals
#' that exceed some ratio * average_residuals some 

filterAnchors <- function(rts, spans, iterFilter, ratio, frac, loess.parameters){
  i = 0
  while(i < iterFilter){
    residuals <- lapply(spans, function(s){
      model <- loess(rty ~ rtx, data = rts, span = s, degree = 1,
                     family = "symmetric", control = loess.parameters)
      
      res = model$residuals
      
      return(res)
    })
    
    residuals = matrix(unlist(residuals), ncol = length(spans)) %>% abs()
    
    thresholds <- colMeans(residuals) * ratio
    
    ###flagging excessively high residuals
    flags <- lapply(1:length(spans), function(i){
      fl <- (residuals[,i] > thresholds[i])
      return(fl)
    }) %>% unlist %>% matrix(ncol = length(spans))
    
    ###deciding which ordered pairs to keep (including min & max)
    fracs = rowSums(flags)/ncol(flags)
    keep = fracs < frac
    keep[c(1,length(keep))] = TRUE
    
    #filtering ordered pairs; need 20 anchors plus both endpoints
    if(sum(keep) > 22 & sum(keep) < length(keep)){          
      rts = rts[keep,]
    } else {
      break
    }
    
    i = i+1
  }
  
  return(rts)
}


#'This is a helper function for the fit_loess function.
#'Performs cross validation to choose the most appropriate span loess for loess fit. 

crossValidationLoess <- function(rts, spans, seed, loess.parameters){
  N = nrow(rts) - 1
  
  #reproducibility seed                                                        
  set.seed(seed)
  
  #creating 10 folds
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
  }) %>% data.frame()
  
  return(cv_errors)
}




##fit_loess
#'
#'
#'
##
fit_loess <- function(object, spans = seq(0.2, 0.4, by = 0.05), ratio = 2, 
                      frac = 0.5, iterFilter = 1, iterLoess = 10, seed = 100)
{
    anchors = object@anchors
  
    if(nrow(anchors) == 0)
        stop("Missing 'anchors' in metabCombiner object. selectAnchors() has not 
              been called yet")
  
    if(nrow(anchors) <= 20)
        stop("Number of anchors too low (less than 20). Examine parameters from 
              selectAnchors() or try using identities.")
  
    ##parameter value checks 
    if(class(ratio)!= "numeric" | ratio < 1){
        warning("Invalid parameter 'ratio'. Must be a number greater than 1. Setting
                value to default(2).")
        ratio = 2
    }
  
    if(class(iterLoess) != "numeric" | round(iterLoess) < 1){
        warning("parameter 'iterLoess' must be at least 1. Setting value to 
                default(10)")
        iterLoess = 10
    }
  
    if(any(spans <= 0 | spans >= 1) | !is.numeric(spans))
        stop("all values in parameter 'spans' must be numeric between 0 & 1")
  
    if(frac <= 0 | frac > 1 | !is.numeric(frac)){
        warning("parameter 'frac' must be numeric between 0 & 1. Setting value to
                 default(0.5)")
        frac = 0.5
    }
  
    loess.parameters = stats::loess.control(iterations = iterLoess, 
                                            surface = "direct")
  
    ##appending minimum and maximum RT values to retention time lists
    
    cTable = combinerTable(object)
    
    minRTx = min(cTable$rtx)
    minRTy = min(cTable$rty)
    maxRTx = max(cTable$rtx)
    maxRTy = max(cTable$rty)
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
  
    object@model = best_model
  
    return(object)
})


