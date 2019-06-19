## Combine unequal names or merge identical strings
#' 
#'
##
combine_strings <- function(s1, s2){
    ifelse(tolower(s1) == tolower(s2), s1, paste(s1,s2, sep = ";"))
}

## Select common pre-identified compounds as anchors.
#'
#'
## 
identityAnchorSelection <- function(cTable, windX, windY, useID){
    if(useID != TRUE)
        return(cTable)
      
    length = nrow(cTable)
    
    ids = combine_strings(cTable$idx, cTable$idy)
    
    mismatchIds = grep(";", ids)
    emptyIds = which(ids == "")
    nonIds = c(emptyIds, mismatchIds)
    
    matches = which(!((1:length) %in% nonIds))
    duplicateIds = matches[duplicated(tolower(cTable$idx[matches]))]
    
    if(length(Ids) > 0){
        cTable$labels[matches] = "I"
        
        cTable$labels[duplicateIds] = "N"
        
        cTable$labels = .Call("selectAnchorsByID",
                              labels = cTable$labels,
                              Ids = Ids,
                              rtx = cTable$rtx,
                              rty = cTable$rty,
                              windX = windX,
                              windY = windY)
    }
    
    return(cTable)
    
}

##iterativeAnchorsSelection
#'
#'
#'
#'
##
iterativeAnchorSelection <- function(cTable, windX, windY, swap = FALSE){
    if(swap)
        cTable$labels = .Call("selectIterativeAnchors",
                               labels = cTable$labels,
                               rtx = cTable$rty,
                               rty = cTable$rtx,
                               Qx = cTable$Qy,
                               Qy = cTable$Qx,
                               windX = windY,
                               windY = windX)
    else
        cTable$labels = .Call("selectIterativeAnchors",
                                labels = cTable$labels,
                                rtx = cTable$rtx,
                                rty = cTable$rty,
                                Qx = cTable$Qx,
                                Qy = cTable$Qy,
                                windX = windX,
                                windY = windY)
      
  
    cTable <- dplyr::filter(cTable, labels != "N")
    
    return(cTable)
}


## selectAnchors()
#'
#'
##
selectAnchors <- function(object, useID = FALSE, tolMZ = 0.005, tolQ = 0.5, 
                          windX = 0.025, windY = 0.025){
    if(class(object) != "metabCombiner")
        stop(base::paste(object, "is not a metabCombiner object"), sep = " ")
  
    if(tolMZ <= 0 | tolQ <= 0 | windX <= 0 | windY <= 0)
        stop("tolMZ, tolQ, windX, windY must all be positive constants")
  
    if(!is.logical(useID))
        stop("Expected logical value for parameter 'useID'")
  
    cTable <- combinerTable(object)
    
    if(nrow(cTable) == 0)
        stop("Empty Combiner table. mzGroup() has not been called yet.")
    
    cTable = dplyr::select(cTable, id, adduct, mzx, mzy,rtx, rty, Qx, Qy, group) %>%
             dplyr::mutate(mzdiff = abs(mzx - mzy), 
                           Qdiff = abs(Qx - Qy), 
                           labels = rep("P", nrow(cTable))) %>%
             identityAnchorSelection(windX = windX, windY = windY, useID = useID) %>% 
             dplyr::filter((abs(mzdiff) < tolMZ & abs(Qdiff) < tolQ) | labels == "I") %>%
             dplyr::filter(labels != "N") %>%
             dplyr::select(-c(mzdiff, Qdiff))
  

    anchorlistXY <- iterativeAnchorSelection(cTable, windX = windX, windY = windY, swap = FALSE)
    anchorlistYX <- iterativeAnchorSelection(cTable, windX = windX, windY = windY, swap = TRUE)
    
    anchorlist <- dplyr::intersect(anchorlistXY, anchorlistYX)
    
    if(nrow(anchorlist) < 20)
        warning("Number of anchors is less than 20. Consider using
                different parameters or using identities.")
    
    object@anchors = anchorlist      
    
    return(object)
}





