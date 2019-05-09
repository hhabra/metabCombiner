##mzdiff
#'
#'
#'
#'
#'
#'
##
mzdiff <- function(rows, mzx, mzy, usePPM){
    diffs = abs(mzx[rows] - mzy[rows])
  
    if(usePPM)
        diffs = 1e6 * diffs / sapply(rows, function(row) min(mzx[row], mzy[row])) 
    
    return(diffs)
}





#'
#'@description used to compare ID or adduct strings
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
#' @param A Numeric Parameter penalizing m/z differences
#' 
#' @param B Numeric Parameter penalizing projected and observed retention times
#' 
#' @param G Numeric Paramater penalizing Q differences
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
#' @return Pairwise score between two grouped features(between 0 & 1) evaluating the 
#'         likelihood of a compound match.
#' 
#' @examples 
#' scorePair(A = 100, B = 5, C = 0.4, mzdiff = 0.003, rtdiff = 0.5, qdiff = 0.1, rtrange = 20, adductMatch = TRUE, adductPenalty = 1)
#' 
scorePairs <- function(A, B, C, mzdiff, rtdiff, qdiff, rtrange, usePPM, 
                       adductMatch, adductPenalty){
   
    nlogScore = A * abs(mzdiff) + B * abs(rtdiff) / rtrange + C * abs(qdiff)
  
    score = base::exp(-nlogScore)
  
    return(score)
}

scorePairs(A = 100, B = 5, C = 0.4, rtrange = 15, mzdiff = abs(cT$mzx - cT$mzy),
           rtdiff = abs(cT$rty - cT$rtProj), qdiff = abs(cT$Qx - cT$Qy),
           adductPenalty = 2)


##calcScores
#'
#'
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
  
    cTab = combinerTable(object)
  
    fit = match.arg(fit)
    model = getModel(object, fit = fit)
    
    rtrange = max(model[["y"]]) - min(model[["y"]])
    
    if(is.null(model))
        stop(paste("object missing model of type ", fit, sep ="" ))
  
    if(useAdduct == FALSE)
        adductPenalty = 1
  
    #handling groups parameter
    if(is.null(groups))
        groups = 1:max(cTab$group)
  
    if(any(!groups %in% cTab$group))
        stop("Invalid parameter 'groups'")
  
    #rows whose scores are computed according to groups
    rows = which(cTab$group %in% groups) 
    
    cTab$rtProj[rows] = predict(model, newdata = cTab$rtx[rows])
    
    cTab$score[rows] = scorePairs(
                        A = A, B = B, C = C,
                        mzdiff = mzdiff(rows, cTab$mzx, cTab$mzy, usePPM),
                        rtdiff = abs(cTab$rty[rows] - cTab$rtProj[rows]),
                        qdiff = abs(cTab$Qx[rows] - cTab$Qy[rows]),
                        rtrange = rtrange,
                        usePPM = usePPM,
                        adductMatch = adductMatch,
                        adductPenalty = adductPenalty
    ) %>% round(digits = 4)

    cTab[rows,] = cTab[rows,] %>% group_by(mzx, rtx) %>% 
                    mutate(rankX_Y = dense_rank(desc(score))) %>%
                    ungroup() %>% group_by(mzy, rty) %>% 
                    mutate(rankY_X= dense_rank(desc(score))) %>%
                    ungroup()
    
    object@combinerTable = arrange(cTab, group, desc(score)) 
    
    return(object)
}


