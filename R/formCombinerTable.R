## Form a combinerTable from computed m/z groups.
#' 
#' @param object
#' 
#' @return 
## 
formCombinerTable <- function(object){
    xdata = getData(object, "x") %>% dplyr::filter(group > 0)
    ydata = getData(object, "y") %>% dplyr::filter(group > 0)
    
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
                                idy = rep(yGrp$id, each = nrow(xGrp)),
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
                                adducty = rep(yGrp$adduct, each = nrow(xGrp)),
                                xsamps, ysamps, 
                                stringsAsFactors = FALSE, check.names = FALSE)
    
        return(groupTable)
    })

    #combine groups into data frame
    object@combinerTable = do.call(rbind.data.frame, combinerTables)
  
    return(object)
}




## Form a combinerTable from computed m/z groups.
#' 
#' @param object
#' 
#' @return 
## 
formCombinerTable <- function(object){
    xdata = getData(object, "x") %>% dplyr::filter(group > 0) %>% arrange(group)
    ydata = getData(object, "y") %>% dplyr::filter(group > 0) %>% arrange(group)
  
    maxGrp = max(xdata[["group"]])
    
    grpCountX = as.integer(table(xdata[["group"]]))
    grpCountY = as.integer(table(ydata[["group"]]))
    totalRows = sum(grpCountX * grpCountY)
    
    xreps = rep(grpCountY, times = grpCountX)

    xCombine = dplyr::slice(xdata, rep(1:n(), times = xreps))
    
    yreps = lapply(1:maxGrp, function(number){
        counts = rep(which(ydata$group == number), times = grpCountX[number])
      
        return(counts)
    }) %>% unlist()
    
    yCombine = ydata[yreps,]
    
    cTable = data.frame(idx = xCombine[["id"]], idy = yCombine[["id"]],
                        mzx = xCombine[["mz"]], mzy = yCombine[["mz"]],
                        rtx = xCombine[["rt"]], rty = yCombine[["rt"]],
                        rtProj = numeric(totalRows),
                        Qx = xCombine[["Q"]], Qy = yCombine[["Q"]],
                        group = xCombine[["group"]],
                        score = numeric(totalRows),
                        rankX_Y = numeric(totalRows),
                        rankY_X = numeric(totalRows),
                        adductx = xCombine[["adduct"]], 
                        adducty = yCombine[["adduct"]],
                        xCombine[,7:ncol(xCombine)],
                        yCombine[,7:ncol(yCombine)],
                        stringsAsFactors = FALSE, check.names = FALSE
    )
                        
  
    #combine groups into data frame
    object@combinerTable = cTable
  
    return(object)
}












