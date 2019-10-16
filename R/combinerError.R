##Contains Error-checking Functions

##
#' @title Obtain Errors for Combiner Object Checks
#'
#' @description This function stores and returns a customized error message
#' when checking the validity of certain objects.
#'
#' @param errNo   integer error code.
#'
#' @param type  character object type (either "combinerTable", "metabCombiner"
#' or "metabData")
#'
#' @return A customized error message for specific object check.
##
combinerError <- function(errNo, type){
    errors = switch(type,
        "combinerTable" =
           c("combinerTable is not a data frame",
              "combinerTable has fewer than expected columns",
              "combinerTable is empty",
              "combinerTable column names different from expected column names",
              "combinerTable column classes different from expected classes",
              "At least one invalid group number in combinerTable"),

        "metabCombiner" =
            c("Input object is not of class 'metabCombiner'",
              "combinerTable in input object is invalid"),

        "metabData" =
            c("Input object is not of class 'metabData'",
              "Empty Data in input object",
              "Data in input object has fewer than expected columns",
              "Data in input object has invalid column names",
              "Data in input object has invalid column classes")
        )

    return(errors[errNo])

}


##
#' @title Determine If Object is a Valid combinerTable
#'
#' @description Checks whether input object is a valid metabData.Returns an
#' integer code if invalid. Function is used alongside \code{combinerError}.
#'
#' @param object  Any R object.
#'
#' @return 0 if object is a valid Combiner Table; an integer code otherwise
##
isCombinerTable <- function(object){
    if(class(object) != "data.frame")
        return(1)

    expected = c("idx", "idy", "mzx", "mzy", "rtx", "rty", "rtProj", "Qx", "Qy",
                 "group", "score", "rankX", "rankY", "adductx", "adducty")

    if(ncol(object) <= length(expected))
        return(2)

    if(nrow(object) == 0)
        return(3)

    if(!identical(names(object)[1:length(expected)], expected))
        return(4)

    coltypes <- as.character(lapply(object, class))
    expected_classes = c("character","character", "numeric","numeric","numeric",
                         "numeric", "numeric","numeric","numeric","integer",
                         "numeric", "integer", "integer","character","character")

    if(!identical(coltypes[1:length(expected)], expected_classes))
        return(5)

    if(any(object[["group"]] <= 0))
        return(6)

    return(0)
}

##
#' @title Determine if object is a valid metabCombiner object
#'
#' @description Checks whether input object is a valid metabCombiner.Returns an
#' integer code if invalid. Function is used alongside \code{combinerError}.
#'
#' @param object Any R object.
#'
#' @return 0 if object is a valid metabData object; an integer code otherwise
##
isMetabCombiner <- function(object){
    if(class(object) != "metabCombiner")
        return(1)

    combiner_table_code = isCombinerTable(combinerTable(object))

    if(combiner_table_code){
        warning(combinerError(combiner_table_code, "combinerTable"))
        return(2)
    }

    return(0)
}



#'
#' @title Determine validity of input metabData object
#'
#' @description Checks whether input object is a valid metabData.Returns an
#' integer code if invalid. Function is used alongside \code{combinerError}.
#'
#' @param object  Any R object
#'
#' @return 0 if object is a valid metabData object; an integer code otherwise.
##
isMetabData <- function(object){
    if(class(object) != "metabData")
        return(1)

    data = getData(object)

    expected = c("id", "mz", "rt", "adduct", "Q", "group")

    if(ncol(data) <= length(expected))
        return(2)

    if(nrow(data) == 0)
        return(3)

    if(!identical(names(data)[1:length(expected)], expected))
        return(4)

    coltypes <- as.character(lapply(data, class))
    expected_classes = c("character","numeric","numeric",
                         "character","numeric", "integer")

    if(!identical(coltypes[1:length(expected)], expected_classes))
        return(5)

    return(0)
}
