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
    
    maxGroup = max(xdata$group)
  
    ##form individual m/z group information
    combinerTables = lapply(1:maxGroup, function(number){
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
    
        groupTable = data.frame(id = id,
                                adduct = adduct,
                                mzx = rep(xgroup$mz, times = nrow(ygroup)),
                                mzy = rep(ygroup$mz, each = nrow(xgroup)),
                                rtx = rep(xgroup$rt, times = nrow(ygroup)),
                                rty = rep(ygroup$rt, each = nrow(xgroup)),
                                rtProj = rep(0, times = nrow(xgroup) * nrow(ygroup)),
                                Qx = rep(xgroup$Q, times = nrow(ygroup)),
                                Qy = rep(ygroup$Q, each = nrow(xgroup)),
                                group = rep(number, times = nrow(xgroup) * nrow(ygroup)),
                                score = rep(0, times = nrow(xgroup) * nrow(ygroup)),
                                rankX_Y = as.vector(rankX_Y),
                                rankY_X = as.vector(rankY_X),
                                xsamps, ysamps, 
                                stringsAsFactors = FALSE, check.names = FALSE)
    
        return(groupTable)
    })

    #combine groups into data frame
    object@combinerTable = do.call(rbind.data.frame, combinerTables)
  
    return(object)
}



