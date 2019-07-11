#'Constructor for the metabData object.
#'
#'@param table      Path to file containing feature table or data.frame object 
#'                  containing features
#'
#'@param mz         Character. Name for column containing m/z values. If "mz",
#'                  will search for {mz, m/z, m.z, mass}
#'
#'@param rt         Character. Name for column containing retention time values. 
#'                  If "rt", will search for {rt, retention time, r.t, RT}.
#'
#'@param id         Character. Name for column containing compound identifiers.
#'                  If "id", will search for {id, compound, feature}
#'
#'@param adduct     Character. Name for column containing adduct labels. If 
#'                  "adduct", will search for {adduct, adducts, annotation}.
#'
#'@param samples    Character. Names of columns containing sample values. If 
#'                  "detect", finds longest stretch of consecutive numeric columns.
#'
#'@param extra      Character. Names of (additional) user-supplied columns.  
#'
#'@param rtmin      Numeric. Minimum retention time for analysis.
#'
#'@param rtmax      Numeric. Maximum retention time for analysis.
#'
#'@param misspc     Numeric. Threshold missingness percentage for analysis.
#'
#'@param zero       Logical. Whether to consider zero values as missing.
#'
#'@param measure    Character. Central abundance measure, either "median" or "mean".
#'
#'@param impute     Logical. Whether to impute value for missing sample measures.
#'
#'@param imputeVal  Imputed value. One of "median", "mean", or numeric value.
#'
#'@param duplicate  Numeric ordered (mz, rt) tolerance pair for duplicate features.   
#'
metabData <- function(table, mz = "mz", rt = "rt", id = "id",
                      adduct = "adduct", samples = "detect", extra = NULL,
                      rtmin = "min", rtmax = "max", misspc = 50, zero = FALSE, 
                      measure = c("median", "mean"), impute = FALSE, imputeVal = 100,
                      duplicate = c(0.003, 0.05))
{  
    if(missing(table))
        stop("required argument 'table' is missing with no default")
    
    if(!is.character(mz))
        stop("non-character argument for variable 'mz'")
    
    if(!is.character(rt))
        stop("non-character argument for variable 'rt'")
    
    if(!is.character(samples))
        stop("non-character argument for variable 'samples'")
    
    if(!is.character(id))
        id = NULL
    
    if(!is.character(adduct))
        adduct = NULL
    
    if(!is.character(extra))
        extra = NULL
    
    if (misspc >= 100 | misspc < 0 | !is.numeric(misspc))
        stop("Parameter 'misspc' must be a numeric value from [0,100)")
    
    if(typeof(table) == "character")
        table <- readData(table)
    else if(class(table) != "data.frame")
        stop("argument 'table' must be a data.frame or path to data file")
    
    if(!is.logical(zero)){
        warning("non-logical value for argument 'zero'. Setting to default (FALSE)")
        zero = FALSE
    }
    
    if(!is.logical(impute)){
        warning("Parameter 'impute' is non-logical. Setting to default (FALSE)")
        impute = FALSE
    }
    
    measure = match.arg(measure)
    
    newData <- methods::new("metabData")
    
    newData = detectFields(Data = newData, table = table,
                           mz = mz, rt = rt, id, adduct = adduct, 
                           samples = samples, extra = extra)
    
    newData = adjustData(Data = newData, misspc = misspc, measure = measure, 
                         rtmin = rtmin, rtmax = rtmax,
                         zero = zero, impute = impute, imputeVal = imputeVal,
                         duplicate = duplicate)
    
    return(newData)
}






