#' @title Read Metabolomics Feature Table From File.
#'
#' @description Reads in the input file as a data frame. Only a tab-delimited
#' .txt file or a .csv file will be accepted.
#'
#' @param file  character. File path to tab-delimited .txt or .csv file
#'
#' @return  The input file read in as a data frame
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

#' @title Detect Keywords in Column Names
#'
#' @description Helper function for metabData() constructor, used to detect
#'   column corresponding to different fields.
#'
#' @param keywords    character regular expression to search for in \code{names}.
#'
#' @param names   character. A vector of column names.
#'
#' @param type    character. One of either "single" to get the first matching
#' column or "multiple" to retrieve all matching columns.
#'
#' @return
#' Integer index or indices of column(s) of keywords; NULL if none found.
#'
#' @noRd
detect <- function(keywords, names, type){
    if(is.null(keywords))
        return(NULL)

    keywords = paste(keywords, collapse = "|")

    indices = grep(keywords, names)

    if(length(indices) == 0)
        return(NULL)

    #searching for a single column, e.g. m/z, rt columns
    if(type == "single")
        return(indices[1])

    #searching for multiple columns, e.g. sample & extra columns
    else
        return(indices)
}

#' @title Automatic Detection of Experimental Sample Columns
#'
#' @description Searches for the longest continuous stretch of numeric columns
#'   to treat as experimental sample columns if no keyword supplied.
#'
#' @param colnames character names of table columns
#'
#' @param coltypes  types of columns
#'
#' @noRd
detectSamples <- function(colnames, coltypes){
    if(all(coltypes == "numeric"))
       return(colnames)

    if (base::all(coltypes != "numeric"))
       stop("no numeric sample columns could be detected")

    consec = base::rle(coltypes)    #looks for consecutive column types

    consec.numeric = (consec[["values"]] == "numeric")
    longest.numeric = base::max(consec$lengths[consec.numeric])

    index = which.max(consec.numeric & consec[["lengths"]] == longest.numeric)

    firstSample = cumsum(consec[["lengths"]])[index] - longest.numeric + 1
    lastSample = cumsum(consec[["lengths"]])[index]

    samples = colnames[firstSample:lastSample]

    return(samples)
}

#' @title Check m/z Column Values
#'
#' @description  Select the column of the input table corresponding to the m/z
#'   values and checks if it matches requirements
#'
#' @param table   An untargeted metabolomics data frame received as input
#'
#' @param col     Column from input data that contains the mz values.
#'
#' @return    m/z value vector
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
#' @description
#' Select column of the input table corresponding to the rt values and perform
#' condition-checks.
#'
#' @param table     An untargeted metabolomics data frame received as input.
#' @param col       Column from input table that contains the rt values
#'
#' @return    numeric retention time value vector
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
#' @description
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
#' @title Process ID & Adduct Columns
#'
#' @description Optional function to process character (ID/ adduct) columns.
#'
#' @param table  An untargeted metabolomics data frame received as input
#'
#' @param col   Column from input data that contains a character field
#'
#' @return   vector of characters
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
#' @description This function ensures that metabolomics datasets used as inputs
#' for the program possess all of the required fields, plus any optional columns
#' that may appear in the final report table.
#'
#' @param Data  a \code{metabData} object.
#'
#' @param table data frame containing metabolomics features or path to
#' metabolomics data file.
#'
#' @param mz    Character name(s) / regular expressions associated with data
#' column containing m/z values. The first column whose name contains this
#' expression will be selected for analysis.
#'
#' @param rt    Character name(s) / regular expression associated with data
#' column containing retention time values. The first column whose name contains
#' this expression will be selected for analysis.
#'
#' @param id     Character name(s) or regular expression associated with data
#' column containing metabolomics feature identifiers. The first column whose
#' name contains this expression will be selected for analysis.
#'
#' @param adduct   Character name(s) or regular expression associated with data
#' column containing adduct, formula, or additional annotations. The first
#' column whose name contains this expression will be selected for analysis.
#'
#' @param samples   Character names of columns containing sample values. All
#' numeric columns containing these keywords are selected for analysis. If no
#' keywords given, will search for longest stretch of numeric columns remaining.
#'
#' @param extra  Character names of columns containing additional feature
#' information, e.g.  non-analyzed sample values. All columns containing these
#' keywords are selected for analysis.
#'
#' @param Q   Character name(s) or regular expression associated with numeric
#' feature abundance quantiles.
#'
#' @return  an initialized and formatted \code{metabData} object.
#'
detectFields <- function(Data, table, mz, rt, id, adduct, samples, extra, Q)
{
    colNames = names(table)

    ##detecting and processing m/z value column
    mzCol = detect(mz, colNames, type = "single")

    if(is.null(mzCol))
        stop("m/z column undefined or used more than once.")

    new_mz <- selectMZ(table = table, col = mzCol)

    table = table[,-mzCol]
    colNames = colNames[-mzCol]

    ##detecting and processing rt value column
    rtCol = detect(rt, colNames, type = "single")

    if(is.null(rtCol))
        stop("rt column undefined or used more than once.")

    new_rt <- selectRT(table, col = rtCol)

    table = table[,-rtCol]
    colNames = colNames[-rtCol]

    ##detecting (optional) identities column; empty column if missing or null
    idCol = detect(id, colNames, type = "single")
    new_id <- selectColumn(table, col = idCol)

    if(!is.null(idCol)){
        table = table[,-idCol]
        colNames = colNames[-idCol]
    }

    ##detecting (optional) adduct column; empty column if missing or null
    adductCol = detect(adduct, colNames, type = "single")
    new_adduct <- selectColumn(table, col = adductCol)

    if(!is.null(adductCol)){
        table = table[,-adductCol]
        colNames = colNames[-adductCol]
    }

    ##detecting (optional) Q column
    QCol = detect(Q, colNames, type = "single")
    new_Q <- selectQ(table, col = QCol)

    if(!is.null(QCol)){
        table = table[,-QCol]
        colNames = colNames[-QCol]
    }

    ##finding (additional, non-analyzed) extra columns
    extraCols = detect(extra, colNames, type = "multiple")

    if(!is.null(extraCols)){
        new_extra = table[,extraCols]
        Data@extra = colNames[extraCols]
        table = table[,-extraCols]
        colNames = colNames[-extraCols]
    }

    else{
        new_extra = NULL
        Data@extra = character(0)
    }

    ##finding analyzed sample columns
    sampleCols = detect(samples, colNames, type = "multiple")

    if(is.null(sampleCols)){
        warning("No column contains keywords in argument 'samples'. Automated
                detection of sample columns may be less accurate.")

        coltypes <- base::as.character(lapply(table, class))

        samples = detectSamples(colNames, coltypes)
    }

    else{
        coltypes <- as.character(lapply(table[,sampleCols], class))

        sampleCols = sampleCols[which(coltypes %in% c("numeric", "integer"))]

        if(length(sampleCols) == 0)
            stop("no numeric columns found using 'samples' argument")

        samples = colNames[sampleCols]
    }

    new_values = table[,samples]

    Data@samples = samples

    Data@data = data.frame(id = new_id, mz = new_mz, rt = new_rt,
                           adduct = new_adduct, Q = new_Q,
                           new_values, check.names = FALSE,
                           stringsAsFactors = FALSE)

    if(length(extra) > 0)
        Data@data[,Data@extra] = new_extra

    return(Data)
}






