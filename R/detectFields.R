#' @title Read Metabolomics Feature Table From File.
#'
#' @description Reads in the input file as a data frame. Only a tab-delimited
#'   .txt file or a .csv file will be accepted.
#'
#' @param file     character. File path to tab-delimited .txt or .csv file
#'
#' @return     The input file read in as a data frame
#'
#' @examples
#' readData("data1.csv")
#' readData("data2.txt")
#'
#' @noRd
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

##
#' @title Detect Keywords in Column Names
#'
#' @description Helper function for metabData() constructor, used to detect
#'   column corresponding to different fields. If one of {"mz", "rt", "id", or
#'   "adduct"}, looks for additional related keywords, or exact name otherwise.
#'
#' @param field    Character. Name of column to look for.
#'
#' @param names    Character vector containing the names of the data
#'
#' @param coltypes    Character vector containing the types for each column
#'
#' @param exclude   Numeric vector of column indices to exclude
#'
#' @return        Integer index of column (if keyword found), or NULL.
#'
#' @noRd
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

#' @title Automatic Detection of Experimental Sample Columns
#'
#' @description Searches for the longest continuous stretch of numeric columns
#'   to treat as experimental sample columns.
#'
#' @param colnames
#'
#' @param coltypes
#'
#' @noRd
detectSamples <- function(colnames, coltypes){
    if(all(coltypes == "numeric"))
       return(colnames)

    if (base::all(coltypes != "numeric"))
       stop("no numeric sample columns could be detected")

    consec = base::rle(coltypes)    #looks for consecutive column types

    consec.numeric = (consec$values == "numeric")
    longest.numeric = base::max(consec$lengths[consec.numeric])

    index = which.max(consec.numeric & consec[["lengths"]] == longest.numeric)

    firstSample = base::cumsum(consec$lengths)[index] - longest.numeric + 1
    lastSample = base::cumsum(consec$lengths)[index]

    samples = colnames[firstSample:lastSample]

    return(samples)
}

#' @title Check m/z Column Values
#'
#' @description  Select the column of the input table corresponding to the m/z
#'   values and checks if it matches requirements
#'
#' @param table     An untargeted metabolomics data frame received as input
#'
#' @param col       Column from input data that contains the mz values.
#'
#' @return          m/z value vector
#'
#' @noRd
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

#' @title Process Retention Time Column
#'
#' Select column of the input table corresponding to the rt values and perform
#' condition-checks.
#'
#' @param table     An untargeted metabolomics data frame received as input.
#' @param col       Column from input table that contains the rt values
#'
#' @return          rt value vector
#'
#' @noRd
selectRT <- function(table, col){
    if(is.null(col))
        stop("retention time column undefined or used more than once.")

    rts <- table[,col]

    if(class(rts) != "numeric")
        stop("retention time column must be numeric")

    if(any(rts <= 0) | any(is.na(rts)))
        stop("retention time values must be positive with no missing values")

    if(any(rts > 100))
        warning("Default program values expect time in minutes.
                Be sure to use correct values if time units are in seconds!")

    return(rts)
}

#' @title Process Quantiles Column
#'
#' Select column of the input table corresponding to the Q values and perform
#' condition-checks. If null, return a vector of 0.
#'
#' @param table     An untargeted metabolomics data frame received as input.
#' @param col       Column from input table that contains the Q values
#'
#' @return          rt value vector
#'
#' @noRd
selectQ <- function(table, col){
    if(is.null(col)){
        Q = rep(0, nrow(table))
        return(Q)
    }

    Q <- table[,col]

    if(!is.numeric(Q))
        stop("Q column must be numeric")

    if(any(is.na(Q)))
        stop("Q column must not have missing values")

    return(Q)
}



##
#' @title
#'
#' @description Optional function to allow for a column of user-supplied IDs.
#'
#' @param table  An untargeted metabolomics data frame received as input
#'
#' @param col   Column from input data that contains a character field
#'
#' @return   vector of ids
#'
#' @examples
#'
#' @noRd
selectColumn <- function(table, col = NULL){
    if (is.null(col)){
        column <-  rep("", nrow(table))
        return(column)
    }

    column <- table[,col] %>% as.character()

    return(column)
}

#' @title Detect metabData Input Columns
#'
#' @param table      Path to file containing feature table or data.frame object
#'                  containing features
#'
#' @param mz         Character. Name for column containing m/z values. If "mz",
#'                  will search for {mz, m/z, m.z, mass}
#'
#' @param rt         Character. Name for column containing retention time values.
#'                  If "rt", will search for {rt, retention time, r.t, RT}.
#'
#' @param id        Character. Name for column containing compound identifiers.
#'                  If "id", will search for {id, compound, feature}
#'
#' @param adduct     Character. Name for column containing adduct labels. If
#'                  "adduct", will search for {adduct, adducts, annotation}.
#'
#' @param samples    Character. Names of columns containing sample values. If
#'                  "detect", finds longest stretch of consecutive numeric columns.
#'
#' @param extra      Character. Names of (additional) user-supplied columns.
#'
#' @return an initialized and formatted \code{metabData} object.
#'
detectFields <- function(Data, table, mz, rt, id, adduct, samples, extra, Q)
{
    #exclude: integer vector ensuring each column is used at most once
    exclude <- integer(5)
    coltypes <- as.character(lapply(table, class))

    ##detecting m/z value column
    mzCol = detect(mz, names = names(table), types = coltypes, exclude = exclude)

    if(is.null(mzCol))
        stop("m/z column undefined or used more than once.")

    new_mz <- selectMZ(table = table, col = mzCol)
    exclude[1] = mzCol        ##updating excluded indices

    ##detecting rt value column
    rtCol = detect(rt, names = names(table), types = coltypes, exclude = exclude)

    new_rt <- selectRT(table, col = rtCol)
    exclude[2] = rtCol        ##updating excluded indices

    ##detecting (optional) identities column; empty column if missing or null
    idCol = detect(id, names = names(table), types = coltypes, exclude = exclude)
    new_id <- selectColumn(table, col = idCol)

    if(!is.null(idCol))
        exclude[3] = idCol

    ##detecting (optional) adduct column; empty column if missing or null
    adductCol = detect(adduct, names = names(table), types = coltypes,
                       exclude = exclude)
    new_adduct <- selectColumn(table, col = adductCol)

    if(!is.null(adductCol))
        exclude[4] = adductCol

    ##detecting (optional) Q column
    QCol = detect(Q, names = names(table), types = coltypes,
                       exclude = exclude)
    new_Q <- selectQ(table, col = QCol)

    if(!is.null(QCol))
        exclude[5] = QCol

    ##removing excluded columns
    table = table[,-exclude]

    ##finding (optional) extra columns
    if(!is.null(extra)){
        if(!any(extra %in% names(table)))
            stop("At least one name in 'extra' field undefined or used
                 more than once")

        new_extra = table[,extra]
        table = table[, !names(table) %in% extra]
    } else
        extra = character(0)

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
                           adduct = new_adduct, Q = new_Q,
                           group = integer(nrow(table)),new_values,
                           check.names = FALSE, stringsAsFactors = FALSE)

    if(length(extra) > 0)
        Data@data[,extra] = new_extra

    return(Data)

}






