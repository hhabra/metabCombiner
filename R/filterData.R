#'filterRT
#'
#' @description  Filters the data in metabData object to a range of retention times 
#'               determined by rtmin & rtmax
#' 
#' @param data   constructed metabolomics data frame.
#'
#' @param rtmin  lower range of retention times for analysis. Defaults to minimum.  
#' 
#' @param rtmax  Upper range of retention times for analysis. Defaults to maximum. 
#'
filterRT <- function(data, rtmin, rtmax){
    rt = data$rt
    
    if(rtmin == "min"){
        rtmin = min(rt)
    }
    
    if(rtmax == "max"){
        rtmax = max(rt)
    }
    
    if(class(rtmin)!= "numeric" | rtmin > rtmax | rtmin < 0 | rtmin > max(rt) |rtmin < min(rt)){
        warning("The supplied 'rtmin value' is invalid. Setting 'rtmin' to minimum observed 
                retention time.")
        rtmin = min(rt)
    }
    
    if(class(rtmax)!= "numeric" | rtmin > rtmax | rtmax < 0 | rtmax > max(rt)){
        warning("The supplied 'rtmax' value is invalid. Setting to maximum observed 
              retention time.")
        rtmax = max(rt)
    }
    
    data = dplyr::filter(data, rt >= rtmin & rt <= rtmax)
    
    return(data)
}

#'imputeVals
#'
#' @description  
#' 
#' @param data   constructed metabolomics data frame.
#'
#' @param samples 
#' 
#' @param imputeVal  

imputeVals <- function(data, samples, imputeVal){
    cols <- which(names(data) %in% samples)
    values = data[,cols]
  
    if (imputeVal == "median"){
        imputedData <- mlr::impute(values, classes = list(
                     numeric = mlr::imputeMedian(),
                     integer = mlr::imputeMedian()))$data
      
        imputedData <- apply(values, 1, function(row){
            imp = median(row, na.rm = TRUE)
            as.numeric(ifelse(is.na(row), imp, row))
        }) %>% t() %>% as.data.frame()
    }else if (imputeVal == "mean"){
        imputedData <- apply(values, 1, function(row){
            imp = mean(row, na.rm = TRUE)
            as.numeric(ifelse(is.na(row), imp, row))
        }) %>% t() %>% as.data.frame()
    }else if(is.numeric(imputeVal) & imputeVal >= 0){
        if(require("mlr"))
            imputedData <- mlr::impute(values, classes = list(
                  numeric = mlr::imputeConstant(imputeVal),
                  integer = mlr::imputeConstant(imputeVal)))$data
        
        else
          imputedData <- apply(values, 1, function(row){
              as.numeric(ifelse(is.na(row), imputeVal, row))
          }) %>% t() %>% as.data.frame()
        
    }else{
        warning("Invalid parameter 'imputeVal'. No imputation performed.")
        
        return(data)
    }  
   
   data[,cols] = imputedData
    
   return(data)
}

#'findDuplicates
#'
#' @description  
#' 
#' @param data   constructed metabolomics data frame.
#'
#' @param counts 
#' 
#' @param missing
#' 
#' @param duplicate 

findDuplicates <- function(data, counts, missing, duplicate){
    if(length(duplicate)!= 2)  
        stop("Parameter 'duplicate' must be a numeric, positive ordered pair (tolMZ, tolRT)")
  
    tolMZ = duplicate[1]
    tolRT = duplicate[2]
  
    if(!is.numeric(tolMZ) | !is.numeric(tolRT))
        stop("Parameter 'duplicate' must be a numeric, positive ordered pair (tolMZ, tolRT)")
    
    if(tolMZ <= 0 | tolRT <= 0)
        stop("Parameter 'duplicate' must be a numeric, positive ordered pair (tolMZ, tolRT)")
    
    dataM = dplyr::select(data, mz,rt) %>% 
            dplyr::mutate(counts = counts, missing = missing, index = 1:nrow(data)) %>%
            dplyr::arrange(mz)                                              
  
    dataM$labels = .Call("findDuplicates", dataM$mz, dataM$rt, nrow(dataM), tolMZ, tolRT, 
                       data$missing, data$counts)
  
    duplicates = dataM$index[dataM$labels == 1]
  
    return(duplicates)
}

#'filterData
#'
#' @description  
#' 
#' @param data   constructed metabolomics data frame.
#'
#' @param samples 
#' 
#' @param misspc
#' 
#' @param measure
#' 
#' @param rtmin
#'   
#' @param rtmax
#' 
#' @param zero
#' 
#' @param impute
#' 
#' @param imputeVal
#' 
#' @param duplicate                                                                         


filterData <- function(data, samples, misspc, measure, rtmin, rtmax, zero,
                       impute, imputeVal, duplicate){
    #parameter checks
    if (misspc >= 100 | misspc < 0 | !is.numeric(misspc))
        stop("Parameter 'misspc' must be a numeric value from [0,100)")
    
    if(!is.logical(zero))
        zero = FALSE
    
    if(!is.logical(impute))
        warning("Parameter 'impute' is non-logical. Setting to default (FALSE)")
  
  
    data = filterRT(data, rtmin = rtmin, rtmax = rtmax)
  
    missingpc <- apply(table, 1, function(row){
        if(zero)
            row[row <= 0] <- NA
        
        sum(is.na(row)) / length(row) * 100
    })
    
    keepIndices = which(missingpc < misspc)
    
    data = data[keepIndices,]
    missingpc = missingpc[keepIndices]
    
    ##imputing missing values; setting all missingpc to 0
    if(impute & any(missingpc > 0)){
        data = imputeVals(data = data, samples = samples, imputeVal = imputeVal)
      
        missingpc = rep(0, nrow(data))
    }
    
    if(measure == "median")
        counts <- apply(data[,samples], 1, median, na.rm = TRUE) %>% as.numeric()    
    else if(measure == "mean")
        counts <- apply(data[,samples], 1, mean, na.rm = TRUE) %>% as.numeric()
    else
        stop("parameter 'measure' must be 'median' or 'mean'")
    
    duplicates = findDuplicates(data = data, counts = counts, missing = missingpc, 
                                duplicate = duplicate)
    
    data = data[-duplicates,]
  
    counts = counts[-duplicates,]
  
  
  
}





