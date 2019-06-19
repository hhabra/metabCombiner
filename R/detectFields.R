#' readData
#'
#' @description Reads in the input file as a data frame. Only a tab-delimited 
#' .txt file or a .csv file will be accepted.
#'
#' @param file     character. File path to tab-delimited .txt or .csv file 
#'  
#' @return         The input file read in as a data frame
#' 
#' @examples       
#' readData("data1.csv")
#' readData("data2.txt")
#' 
readData <- function(file){
    type = substr(file , start = (nchar(file)-3), stop = nchar(file))
    
    if(type == ".csv") 
          data <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE) 
    else if (type == ".txt")
          data <- read.delim(file, sep = "\t", stringsAsFactors = FALSE, 
                             check.names = FALSE) 
    else
         stop("Invalid File Type. Must be a .csv or .txt file.")
    
    
    return(data)
}


#' detect
#'
#' @description Helper function for metabData() constructor, used to detect 
#' column corresponding to different fields. If one of {"mz", "rt", "id", or
#' "adduct"}, looks for additional related keywords, or exact name otherwise. 
#'
#' @param field    Character. Name of column to look for.
#' 
#' @param names    Character vector containing the names of the data
#' 
#' @param coltypes    Character vector containing the types for each column
#' 
#' @param exclude   Numeric vector of column indices to exclude
#' 
#' @return        Name for column pertaining to field, or NULL
#' 
detect <- function(field, names, types, exclude){
    if(is.null(field))
        return(NULL)
    
    snames = character(0)     
    field = as.character(field[1])
    keywords = base::switch(field,  "mz" = c("mz", "m/z", "m.z","MZ", "mass"),
                                    "rt" = c("rt","retention time", "r.t", "RT"),
                                    "id" = c("id", "ID", "compound", "Compound", 
                                             "feature", "Feature"),
                                  "adduct" = c("adduct", "adducts", "annotation",
                                               "Adduct", "Adducts")
    )
    
    N = length(names)
    
    #limit to numeric and non-excluded indices
    if(!is.null(keywords)){
        if(field %in% c("mz","rt"))
            snames = names[types == "numeric" & !(1:N %in% exclude)]
        
        else if(field %in% c("id","adduct")){
            snames = names[types == "character" & !(1:N %in% exclude)]
        }
    }else{  #detect user-given names
        keywords = field
        snames = names[!(1:N %in% exclude)]
    }
    
    #collect first name that matches keyword(s)
    if (any(snames %in% keywords)){
        vname = snames[which(snames %in% keywords)[1]]
        index = which(names == vname)
        return(index)
    }else
        return(NULL)
}

#'detectSamples
#'
#' @param colnames
#' 
#' @param coltypes
#' 
#' @return 
#' 
detectSamples <- function(colnames, coltypes){
    if(all(coltypes == "numeric"))
       return(colnames)
  
    if (base::all(coltypes != "numeric"))
       stop("no numeric sample columns could be detected")
  
    consec = base::rle(coltypes)    #looks for consecutive column types
  
    consec.numeric = (consec$values == "numeric")
    longest.numeric = base::max(consec$lengths[consec.numeric])
    
    index = base::which.max(consec.numeric & consec$lengths == longest.numeric)
    
    firstSample = base::cumsum(consec$lengths)[index] - longest.numeric + 1
    lastSample = base::cumsum(consec$lengths)[index]
    
    samples = colnames[firstSample:lastSample]
    
    return(samples)
}

#' selectMZ
#' 
#' @description  Select the column of the input table corresponding to the m/z 
#' values and checks if it matches requirements
#'
#' @param table     An untargeted metabolomics data frame received as input
#' 
#' @param col       Column from input data that contains the mz values
#' 
#' @return          m/z value vector
#' 
#' @examples
selectMZ <- function(table, col){
    mzs <- table[,col]
    
    if(class(mzs) != "numeric") 
        stop("m/z column must be numeric")
    
    if(any(mzs <= 0)) 
        stop("all m/z values must be strictly positive")
    
    if(any(mzs >= 2000) | any(mzs <= 50))
        warning("Unusual m/z values detected. Be sure that m/z column
                 is correctly chosen.")
    
    return(mzs)
}

#' selectRT
#' 
#' Select column of the input table corresponding to the rt values and perform 
#' condition-checks. Return rt column
#'
#' @param table     An untargeted metabolomics data frame received as input.
#' @param col       Column from input table that contains the mz values
#' 
#' @return          rt value vector
#' 
#' @examples

selectRT <- function(table, col){
    rts <- table[,col]
    
    if(class(rts) != "numeric")  
        stop("retention time column must be numeric")
    
    if(any(rts <= 0)) 
        stop("retention time values must be positive")
    
    if(any(rts > 100))
        warning("Default retention time program values are in minutes. Be sure to use correct values
                if time units are in seconds!")
    
    return(rts)
}  


#' Optional function to allow for a column of user-supplied IDs.
#'
#' @param table      An untargeted metabolomics data frame received as input
#'
#' @param col       Column from input data that contains the IDs
#' 
#' @return          vector of ids
#' 
#' @examples
selectID <- function(table, col = NULL){
    if (!is.numeric(col)){
        ids <-  rep("", nrow(table))
        return(ids)
    }
    
    ids <- table[,col] %>% as.character()
    
    return(ids)
}


#' Optional function to allow for a column of user-supplied annotations.
#'
#' @param table     An untargeted metabolomics data frame received as input
#' @param col       Column from input table that contains the IDs
#' 
#' @return          vector of adducts
#' 
selectAdduct <- function(table, col = NULL){
    if (!is.numeric(col)){
        adducts <-  rep("", nrow(table))
        return(adducts)
    }
    
    adducts <- table[,col] %>% as.character()
    
    return(adducts)
}


#'metabData
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
#'  
#'  
detectFields <- function(Data, table, mz, rt, id, adduct, samples, extra)
{  
    #exclude: integer vector ensures that each column is used at most one time
    exclude <- integer(ncol(table))
    coltypes <- as.character(lapply(table, class))
    
    ##detecting m/z value column
    mzCol = detect(mz, names = names(table), types = coltypes, exclude = exclude)
    
    if(is.null(mzCol))
        stop("m/z column undefined or used more than once.")
    
    new_mz <- selectMZ(table = table, col = mzCol)
    exclude[1] = mzCol        ##updating excluded indices
    
    ##detecting rt value column
    rtCol = detect(rt, names = names(table), types = coltypes, exclude = exclude)
    
    if(is.null(rtCol))
        stop("Retention Time column undefined or used more than once.")
    
    new_rt <- selectRT(table, col = rtCol)
    exclude[2] = rtCol        ##updating excluded indices
    
    ##detecting (optional) identities column; empty column if missing or null
    idCol = detect(id, names = names(table), types = coltypes, exclude = exclude)
    new_id <- selectID(table, col = idCol)
    
    if(!is.null(idCol))
        exclude[3] = idCol
    
    ##detecting (optional) adduct column; empty column if missing or null
    adductCol = detect(adduct, names = names(table), types = coltypes,
                       exclude = exclude)
    
    new_adduct <- selectAdduct(table, col = adductCol)
    
    if(!is.null(adductCol))
        exclude[4] = adductCol
    
    ##removing excluded columns
    table = table[,-exclude]
    
    ##finding (optional) extra columns
    if(is.null(extra)){
        new_extra = base::rep("", nrow(table))
        extra = character(0)
    }else{
        if(!any(extra %in% names(table)))
            stop("At least one name in 'extra' field undefined or used more than once")
      
        new_extra = table[,extra]
        table = table[, !names(table) %in% extra]
    }
    
    Data@extra = extra
    
    ##finding sample value columns
    if(base::identical(samples, "detect")){
        coltypes <- base::as.character(lapply(table, class))
      
        samples = detectSamples(names(table), coltypes)
    }else{
        if(!any(samples %in% names(table)) | any(duplicated(samples)))
            stop("At least one name in 'samples' undefined or non-unique")
    }
    
    new_values = table[,samples]
    
    apply(new_values, 2, function(x){
        if(!(class(x) %in% c("integer", "numeric")))
           stop("at least one selected sample column is non-numeric")
        
    })
    
    Data@samples = samples
    
    Data@data = data.frame(id = new_id, mz = new_mz, rt = new_rt, 
                           adduct = new_adduct, Q = rep(0, nrow(table)), 
                           group = rep(0, nrow(table)),new_values, 
                           check.names = FALSE, stringsAsFactors = FALSE)
    
    if(length(extra) > 0)
        newData@data[,extra] = new_extra
    
    return(newData)
}





