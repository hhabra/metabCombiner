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
    type <- substr(file , start = (nchar(file)-3), stop = nchar(file))

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
#' @param keywords character regular expression to search for in \code{names}.
#'
#' @param names  character. A vector of column names.
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
    keywords <- paste(keywords, collapse = "|")
    indices <- grep(keywords, names)
    if(length(indices) == 0)
        return(NULL)
    if(type == "single")
        return(indices[1])
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
    coltypes <- base::gsub("integer", "numeric", coltypes)
    if(all(coltypes == "numeric"))
        return(colnames)
    if (base::all(coltypes != "numeric"))
        stop("no numeric sample columns could be detected")
    consec <- base::rle(coltypes)    #looks for consecutive column types
    consec.numeric <- (consec[["values"]] == "numeric")
    longest.numeric <- base::max(consec$lengths[consec.numeric])

    index <- which.max(consec.numeric & consec[["lengths"]] == longest.numeric)
    firstSample <- cumsum(consec[["lengths"]])[index] - longest.numeric + 1
    lastSample <- cumsum(consec[["lengths"]])[index]

    samples <- colnames[seq(firstSample, lastSample)]
    return(samples)
}

#' @title Check m/z Column Values
#'
#' @description  Select the column of the input table corresponding to the m/z
#'   values and checks if it matches requirements
#'
#' @param table An untargeted metabolomics data frame received as input
#'
#' @param col Column from input data that contains the mz values.
#'
#' @return  m/z value vector
#'
#' @noRd
selectMZ <- function(table, col){
    if(is.null(col))
        stop("m/z column undefined or used more than once")

    mzs <- table[,col]

    if(!is.numeric(mzs))
        stop("m/z column must be numeric")

    if(any(mzs <= 0) | any(is.na(mzs)))
        stop("all m/z values must be positive and not missing")

    if(any(mzs >= 2000) | any(mzs <= 20))
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
#' @param table An untargeted metabolomics data frame received as input.
#' @param col  Column from input table that contains the rt values
#'
#' @return    numeric retention time value vector
#'
#' @noRd
selectRT <- function(table, col){
    if(is.null(col))
        stop("retention time column undefined or used more than once")

    rts <- table[,col]

    if(!is.numeric(rts))
        stop("retention time column must be numeric")

    if(any(rts <= 0) | any(is.na(rts)))
        stop("retention time values must be positive and not missing")

    if(any(rts > 500))
        warning("program defaults expect time in minutes; ",
                "adjust values if time units are in seconds")
    return(rts)
}

#' @title Process Quantiles Column
#'
#' @description
#' Select column of the input table corresponding to the Q values and perform
#' condition-checks. If null, return a vector of 0.
#'
#' @param table An untargeted metabolomics data frame received as input.
#' @param col   Column from input table that contains the Q values
#'
#' @return  rt value vector
#'
#' @noRd
selectQ <- function(table, col){
    if(is.null(col)){
        Q <- rep(0, nrow(table))
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
#' @return  vector of characters
#'
#' @noRd
#'
#'
selectColumn <- function(table, col = NULL){
    if (is.null(col)){
        column <- rep("", nrow(table))
        return(column)
    }
    column <- as.character(table[,col])
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
#' @param rt Character name(s) / regular expression associated with data
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
#' keywords given, searches for longest stretch of numeric columns remaining.
#'
#' @param extra  Character names of columns containing additional feature
#' information, e.g.  non-analyzed sample values. All columns containing these
#' keywords are selected for analysis.
#'
#' @param Q  Character name(s) or regular expression associated with numeric
#' feature abundance quantiles.
#'
#' @return  an initialized and formatted \code{metabData} object.
#'
detectFields <- function(Data, table, mz, rt, id, adduct, samples, extra, Q)
{
    mzCol <- detect(mz, names(table), type = "single")
    new_mz <- round(selectMZ(table = table, col = mzCol),4)
    table <- table[-mzCol]

    rtCol <- detect(rt, names(table), type = "single")
    new_rt <- round(selectRT(table, col = rtCol),4)
    table <- table[-rtCol]

    idCol <- detect(id, names(table), type = "single")
    new_id <- selectColumn(table, col = idCol)
    if(!is.null(idCol)) table <- table[-idCol]

    adductCol <- detect(adduct, names(table), type = "single")
    new_adduct <- selectColumn(table, col = adductCol)
    if(!is.null(adductCol)) table <- table[-adductCol]

    QCol <- detect(Q, names(table), type = "single")
    new_Q <- selectQ(table, col = QCol)
    if(!is.null(QCol)) table <- table[-QCol]

    sampleCols <- detect(samples, names(table), type = "multiple")
    if(is.null(sampleCols)){
        warning("no column(s) found for argument 'samples'; ",
                "samples column detection may be inaccurate.")
        coltypes <- base::as.character(vapply(table, class, character(1)))
        samples <- detectSamples(names(table), coltypes)
    } else{
        coltypes <- vapply(table, class, character(1))[sampleCols]
        sampleCols <- sampleCols[which(coltypes %in% c("numeric", "integer"))]
        if(length(sampleCols) == 0)
            stop("no numeric columns found using 'samples' argument")
        samples <- names(table)[sampleCols]
    }
    data <- data.frame(rowID = seq_along(new_mz), id = new_id, mz = new_mz,
                        rt = new_rt, adduct = new_adduct, Q = new_Q,
                        table[samples], check.names = FALSE,
                        stringsAsFactors = FALSE)
    table <- dplyr::select(table, -dplyr::all_of(samples))
    extraCols <- detect(extra, names(table), type = "multiple")
    extra <- names(table)[extraCols]
    if(!is.null(extraCols))
        data[extra] <- table[extraCols]
    Data <- update_md(Data, data = data, samples = samples, extra = extra)
    return(Data)
}






