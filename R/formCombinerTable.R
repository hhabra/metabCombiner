## combine_strings
#' 
#' For combining unequal strings or merging identical strings
#'
#'
##
combine_strings <- function(s1, s2){
    ifelse(s1 == s2, s1, paste(s1,s2, sep = ";"))
}

##formCombinerTable
#' 
#' @param object
#' 
#' @return 
## 
formCombinerTable <- function(object){
    xdata = getData(object, "x") %>% filter(group > 0)
    ydata = getData(object, "y") %>% filter(group > 0)
    
    maxGrp = max(xdata$group)
  
    ##form individual m/z group information
    combinerTables = lapply(1:maxGrp, function(number){
        xGrp = dplyr::filter(xdata, group == number)
        yGrp = dplyr::filter(ydata, group == number)
    
        rankX_Y = rep(1, times = nrow(xGrp) * nrow(yGrp))
        rankY_X = rep(1, times = nrow(xGrp) * nrow(yGrp))
    
        xsamps = xGrp[,7:ncol(xGrp)] %>% slice(rep(1:n(), times = nrow(yGrp))) %>% 
                 as.data.frame()
        ysamps = yGrp[,7:ncol(yGrp)] %>% slice(rep(1:n(), each = nrow(xGrp))) %>% 
                 as.data.frame()
    
        groupTable = data.frame(idx = rep(xGrp$id, times = nrow(yGrp)),
                                idy = rep(yGrp$id, each = nrow(yGrp)),
                                mzx = rep(xGrp$mz, times = nrow(yGrp)),
                                mzy = rep(yGrp$mz, each = nrow(xGrp)),
                                rtx = rep(xGrp$rt, times = nrow(yGrp)),
                                rty = rep(yGrp$rt, each = nrow(xGrp)),
                                rtProj = rep(0, times = nrow(xGrp)*nrow(yGrp)),
                                Qx = rep(xGrp$Q, times = nrow(yGrp)),
                                Qy = rep(yGrp$Q, each = nrow(xGrp)),
                                group = rep(number, times = nrow(xGrp)*nrow(yGrp)),
                                score = rep(0, times = nrow(xGrp)*nrow(yGrp)),
                                rankX_Y = as.vector(rankX_Y),
                                rankY_X = as.vector(rankY_X),
                                adductx = rep(xGrp$adduct, times = nrow(yGrp)),
                                adducty = rep(xGrp$adduct, each = nrow(yGrp)),
                                xsamps, ysamps, 
                                stringsAsFactors = FALSE, check.names = FALSE)
    
        return(groupTable)
    })

    #combine groups into data frame
    object@combinerTable = do.call(rbind.data.frame, combinerTables)
  
    return(object)
}



