## compare_strings
#'
## combine_strings
#' 
#' for combining names together
#'
##
combine_strings <- function(s1, s2){
    ifelse(s1 == s2, s1, paste(s1,s2, sep = ";"))
}

##formCombinerTable
#' 
#' @param object
#' 
#' @param maxGroup
#' 
#' @return 
#' 
formCombinerTable <- function(object, maxGroup){
    xdata = getData(object, "x")
    ydata = getData(object, "y")
    
    maxGroup = max(xdata$group)
    
    combinerTables = lapply(1:maxGroup, function(number)){
        xgroup = dplyr::filter(xdata, group == number)
        ygroup = dplyr::filter(ydata, group == number)
        
        rankX_Y = rep(1, times = nrow(xgroup) * nrow(ygroup))
        rankY_X = rep(1, times = nrow(xgroup) * nrow(ygroup))
    
        xsamps = xgroup[,7:ncol(xgroup)] %>% slice(rep(1:n(), times = nrow(ygroup))) %>% 
                 as.data.frame()
        ysamps = ygroup[,7:ncol(ygroup)] %>% slice(rep(1:n(), each = nrow(xgroup))) %>% 
                 as.data.frame()
    
        id = as.vector(outer(xgroup$id, ygroup$id, combine_strings))
        adduct = as.vector(outer(xgroup$adduct, ygroup$adduct, combine_strings))
        
        groupCombinerTable = data.frame(
                                    id = id,
                                    adduct = adduct,
                                    mzx = rep(xgroup$mz, times = nrow(ygroup)),
                                    mzy = rep(ygroup$mz, each = nrow(xgroup)),
                                    rtx = rep(xgroup$rt, times = nrow(ygroup)),
                                    rty = rep(ygroup$rt, each = nrow(xgroup)),
                                    rtProj = rep(0, times = nrow(ygroup)),
                                    Qx = rep(xgroup$Q, times = nrow(ygroup)),
                                    Qy = rep(ygroup$Q, each = nrow(xgroup)),
                                    group = rep(number, times = nrow(xgroup) * nrow(ygroup)),
                                    score = as.vector(scoreMatrix),
                                    rankX_Y = as.vector(rankX_Y),
                                    rankY_X = as.vector(rankY_X),
                                    xsamps, ysamps, 
                                    stringsAsFactors = FALSE, check.names = FALSE)
    
        return(groupCombinerTable)
    })

    object@combinerTable = do.call(rbind.data.frame, combinerTables)
  
    return(object)
}


