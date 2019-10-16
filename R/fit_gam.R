## 
#' @title Filter Outlier Ordered Pairs
#' 
#' @description Helper function for fit_gam.
#' 
#' @param rts 
#' 
#' @param k
#' 
#' @param bs
#' 
#' @param m
#' 
#' @param family
#' 
#' @param method
#' 
#' @param iterFilter
#' 
#' @param optim
#' 
#' @param ratio
#' 
#' @param frac
#'
##
filterAnchorsGAM <- function(rts, k, bs, m, family, method, iterFilter, optim, 
                             ratio, frac,...){
    iteration = 0
  
    while(iteration < iterFilter){
        cat("Performing filtering iteration: ", iteration + 1, "\n", sep = "")
        residuals <- sapply(k, function(ki){
              model <- mgcv::gam(rty ~ s(rtx, k = ki, bs = bs, m = m,...), 
                             data = rts, family = family, method = method, 
                             weights = rts$weights, 
                             optimizer = c("outer", optim), ...)
      
              res = abs(model$fitted.values - model$y)
              return(res)
        })
        
        include = which(rts$weights == 1) 
    
        ##flagging excessively high residuals
        thresholds <- ratio * colMeans(residuals[include,] %>% as.matrix())
    
        flags <- sapply(1:length(k), function(i){
            fl <- (residuals[,i] > thresholds[i])
            return(fl)
        })
    
        ##remove features with high proportion of excessive residuals
        fracs = rowSums(flags)/ncol(flags)
        remove = fracs > frac
        remove[rts$labels == "I"] = FALSE
    
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
#' @description Helper function for \code{fit_gam()}. Determines optimal value
#' of \code{k} basis functions (among user-defined choices) for Generalized
#' Additive Model, using a 10-fold cross validation. 
#'
#' @param rts
#' 
#' @param k
#'
#' @param bs
#' 
#' @param m
#' 
#' @param family
#' 
#' @param method
#' 
#' @param optim
#' 
#' @param seed integer. Psuedo-random seed generator for cross-validation.
#'
#' @return Optimal value for \code{k} as determined by 10-fold cross validation.
#' 
##
crossValidationGAM <- function(rts, k, bs, m, family, method, optim, seed,...){
    
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
                           optimizer = c("outer",optim), 
                           weights = rts_train$weights,
                           ...)
      
            preds <- predict(model, newdata = rts_test)
      
            MSE = sum((preds - rts_test$rty)^2)/ length(preds)
      
            return(MSE)
        })
    
        return(errors)
    })
    
    mean_cv_errors = rowMeans(matrix(cv_errors, nrow = length(k)))
  
    best_k = k[which.min(mean_cv_errors)]
    
    return(best_k)
}

##
#' @title Fit GAM model between retention times
#' 
#' @description Fit a Generalized Additive Model between selected ordered retention 
#' time pairs
#'
#' @param object  a metabCombiner object.
#'
#' @param useID  logical. Option to use matched IDs to inform modeling process.
#'
#' @param k   integer vector. Selection of values to control the number of basis
#' functions for gam construction. Best value chosen by 10-fold cross validation.
#' 
#' @param ratio numeric. Residual multiplier for determining outliers.
#' 
#' @param iterFilter integer. Number of residual filtering iterations to perform.
#' 
#' @param bs   character. Choice of spline method from mgcv; see: ?smooth.terms
#'
#' @param family character. Choice of mgcv family; see: ?family.mgcv
#' 
#' @param m   integer. Basis and penalty order for GAM; see ?mgcv::s
#' 
#' @param method  character. Smoothing parameter estimation method; see: ?gam
#'
#' @param optimizer character. Method to optimize smoothing parameter; see: ?gam
#'
#' @param seed integer. Psuedo-random seed generator for cross validation. 
#' 
#' @return Updated metabCombiner object with a GAM model
## 
fit_gam <- function(object, useID = FALSE, k = seq(10,25, by = 5), iterFilter = 2,
                    ratio = 2, frac = 0.5, bs = c("bs", "ps", "tp"),
                    family = c("scat", "gaussian"), weights = 1, m = c(3,2),
                    method = "REML", optim = "newton", seed = 100, ...)
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
    
    if(any(k <= 2 | k >= nrow(anchors)) | !is.numeric(k))
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
    
    if(!useID)
        anchors$labels = "A"
  
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
                                            optim = optim,
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
                                                  optim = optim, 
                                                  seed = seed,...))
                                    
  
    
    best_model <- mgcv::gam(rty ~ s(rtx, k = best_k, bs = bs, m = m, ...),
                                  data = rts, family = family, method = method, 
                                  optimizer = c("outer", optim), 
                                  weights = rts[["weights"]], ...)
    
    object@model[["gam"]] = best_model
    object@stats[["best_k"]] = best_k
    
    return(object)
} 

