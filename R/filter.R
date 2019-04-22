#'filterRT
#'
#' @description  Filters the data in a metabData object to a range of retention times 
#'               determined by rtmin & rtmax
#' 
#' @param data  metabData object.
#'
#' @param rtmin  lower range of retention times for analysis. Defaults to minimum.  
#' 
#' @param rtmax  Upper range of retention times 
#'

filterRT <- function(data, rtmin, rtmax){
    rt = data$rt
    
    if(rtmin == "min"){
      rtmin = min(rt)
    }
    
    if(rtmax == "max"){
      rtmax = max(rt)
    }
    
    if(class(rtmin)!= "numeric" | rtmin > rtmax | rtmin < 0 | rtmin > max(rt) |
       rtmin < min(rt)){
      warning("The supplied rtmin is invalid. Setting 'rtmin' to minimum observed 
              retention time.")
      rtmin = min(data$rt)
    }
    
    if(class(rtmax)!= "numeric" | rtmin > rtmax | rtmax < 0 | rtmax > max(rt)){
      warning("The rtmax retention time is invalid. Setting to maximum observed 
              retention time.")
      rtmax = max(rt)
    }
    
    object@data = dplyr::filter(object@data, rt >= rtmin & rt <= rtmax)
    
    return(object)
})








filterData <- function(){
  
  
  
  
  
  
  
  
  
}





