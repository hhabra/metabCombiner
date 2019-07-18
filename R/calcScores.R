## calculate differences in m/z between complementary features
#'
#' @param rows integer vector. Indices to calculate m/z differences.
#'
#' @param mzx  numeric vector. m/z values from dataset X.
#'
#' @param mzy  numeric vector. m/z values from dataset Y.
#' 
#' @param usePPM logical. Option to calculate relative differences (ppm).
#' 
#' @return Vector of numeric relative or absolute m/z differences.
##
mzdiff <- function(rows, mzx, mzy, usePPM){
    diffs = abs(mzx[rows] - mzy[rows])
  
    if(usePPM)
        diffs = 1e6 * diffs/ sapply(rows, function(r) min(mzx[r], mzy[r])) 
    
    return(diffs)
}

##
#'
#'
#'
#'
##
compare_strings <- function(string){
    string = base::strsplit("deekee", split = ";")
  
    if(s1 == "" | s2 == "")
       return(TRUE)
  
    if(s1 != s2)
        return(FALSE)
  
  
     return(TRUE)
}

#' Calculate the pairwise score between two grouped features
#' 
#' @param A Numeric. Weight for m/z differences.
#' 
#' @param B Numeric. Weight for differences btw. fitted & observed 
#' retention times
#' 
#' @param G Numeric. Weight for differences in Q (relative abundance).
#'   
#' @param mzdiff Difference between feature m/z values; can be absolute or
#'               relative (ppm)
#' 
#' @param rtdiff Difference between model-projected retention time value ('xdata') 
#'               & observed retention time ('ydata')
#' 
#' @param qdiff Difference between feature quantile Q values
#' 
#' @param adductPenalty Numeric. If useAdduct set to True, used to divide total score
#'                      if adduct labels are mismatched.
#' 
#' @return Pairwise score between two grouped features(between 0 & 1) evaluating
#'         the likelihood of a compound match.
#' 
#' @examples 
#' scorePair(A = 100, B = 5, C = 0.4, mzdiff = 0.003, rtdiff = 0.5, qdiff = 0.1,
#' rtrange = 20, adductMatch = FALSE, adductPenalty = 1)
#' 
scorePairs <- function(A, B, C, mzdiff, rtdiff, qdiff, rtrange, usePPM, 
                       adductMatch, adductPenalty){
   
    nlogScore = A * abs(mzdiff) + B * abs(rtdiff) / rtrange + C * abs(qdiff)
  
    score = base::exp(-nlogScore)
  
    return(score)
}


## Assign pairwise confidence scores between all grouped metabolomics features.
#'
#' @param object  metabCombiner object.
#' 
#' @param A Numeric. Weight for m/z differences.
#' 
#' @param B Numeric. Weight for differences btw. fitted & observed 
#' retention times
#' 
#' @param G Numeric. Weight for differences in Q (relative abundance).
#'
#'
##
calcScores <- function(object, A, B, C, fit = c("loess", "gam"), usePPM = FALSE, 
                       useAdduct = FALSE, adductPenalty = 2, groups = NULL)
{
    if(class(object) != "metabCombiner")
        stop(base::paste(object, "is not a metabCombiner object"), sep = " ")
      
    if(class(A) != "numeric" | class(B) != "numeric" | class(C) != "numeric")
        stop("parameters 'A', 'B', 'C' must be numeric constants")
  
    if(A < 0 | B < 0 | C < 0)
        stop("parameters 'A', 'B', 'C' must be non-negative")
  
    cTable = combinerTable(object)
    
    if(nrow(cTable) == 0)
        stop("missing Combiner table")
    
    fit = match.arg(fit)
    model = getModel(object, fit = fit)
    
    rtrange = max(model[["y"]]) - min(model[["y"]])
    
    if(is.null(model))
        stop(paste("object missing model of type ", fit, sep ="" ))
  
    if(useAdduct == FALSE)
        adductPenalty = 1
  
    #handling groups parameter
    if(is.null(groups))
        groups = 1:max(cTable["group"])
  
    if(any(!groups %in% cTable$group))
        stop("Invalid parameter 'groups'")
  
    #rows whose scores are computed according to groups
    rows = which(cTable$group %in% groups) 
    
    cTable$rtProj[rows] = predict(model, newdata = cTable[rows,])
    
    cTable$score[rows] = scorePairs(
                         A = A, B = B, C = C,
                         mzdiff = mzdiff(rows, cTable$mzx, cTable$mzy, usePPM),
                         rtdiff = abs(cTable$rty[rows] - cTable$rtProj[rows]),
                         qdiff = abs(cTable$Qx[rows] - cTable$Qy[rows]),
                         rtrange = rtrange,
                         usePPM = usePPM,
                         adductMatch = adductMatch,
                         adductPenalty = adductPenalty
    )
    
    cTable[c("rtProj", "score")] = round(cTable[c("rtProj", "score")], 4)
    
    
    ##bug fix for duplicated column names
    cT = cTable[ , !duplicated(colnames(cTable))]

    ##calculate score ranking of features 
    cT[rows,] = cT[rows,] %>% group_by(mzx, rtx) %>% 
                  mutate(rankX_Y = dense_rank(desc(score))) %>%
                  ungroup() %>% group_by(mzy, rty) %>% 
                  mutate(rankY_X= dense_rank(desc(score))) %>%
                  ungroup()
    
    cTable[c("rankX_Y", "rankY_X")] = cT[c("rankX_Y", "rankY_X")]

    object@combinerTable = cTable[with(cTable, order(group, desc(score))), ]
    
    return(object)
}


