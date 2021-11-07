
#' Find Lost Features
#'
#'
#'
#' @noRd
find_lost_features <- function(fdata, forig, type, di)
{
    mzdat <- fdata[grep("^mz", names(fdata), value = TRUE)[di]]
    rtdat <- fdata[grep("^rt", names(fdata), value = TRUE)[di]]
    fdata <- round(cbind.data.frame(mzdat, rtdat),4)
    forig <- round(forig[grep("^mz|^rt", names(forig), value = TRUE)], 4)
    names(forig) <- names(fdata)
    row.names(forig) <- as.character(seq(1,nrow(forig)))
    featdiff <- setdiff(forig, fdata)
    lostfeats <- which(row.names(forig) %in% row.names(featdiff))
    return(lostfeats)
}

#' Fill Feature Data
#'
#' @noRd
fill_featdata <- function(fdata, forig, lostfeats, type, di)
{
    forig <- forig[lostfeats,]
    lost <- seq(nrow(fdata) - length(lostfeats) + 1, nrow(fdata))
    start <- length(grep(paste("^", type, ""), fdata$rowID))

    fdata$rowID[lost] <-  sprintf("%s.0.%s", type,
                                  seq(start+1,start+length(lostfeats)))
    idi <- grep("^id", names(fdata))[di]
    mzi <- grep("^mz", names(fdata))[di]
    rti <- grep("^rt", names(fdata))[di]
    Qi <- grep("^Q", names(fdata))[di]
    adducti <- grep("^adduct", names(fdata))[di]

    fdata[lost, idi] <- forig[grep("^id", names(forig), value = TRUE)]
    fdata[lost, mzi] <- forig[grep("^mz", names(forig), value = TRUE)]
    fdata[lost, rti] <- forig[grep("^rt", names(forig), value = TRUE)]
    fdata[lost, Qi] <- forig[grep("^Q", names(forig), value = TRUE)]
    fdata[lost, adducti] <- forig[grep("^adduct",names(forig),value = TRUE)]
    return(fdata)
}


#' Fill Combined Table
#'
#' @noRd
fill_combinedTable <- function(cTable, origTable, forig, lostfeats,
                               datasets, type, object)
{
    if(type == "x") di <- which(datasets$x == x(object))
    else if(type == "y") di <- which(datasets$y == y(object))
    origTable <- origTable[lostfeats,]
    forig <- forig[lostfeats,]
    lost <- seq(nrow(cTable) - length(lostfeats) + 1, nrow(cTable))

    start <- length(grep(paste("^", type, ""), cTable$rowID))
    cTable$rowID[lost] <-  sprintf("%s.0.%s", type,
                                  seq(start+1,start+length(lostfeats)))
    fields <- paste(c("id", "mz", "rt", "Q", "adduct"), type, sep = "")
    idi <- grep("^id", names(forig))[di]
    mzi <- grep("^mz", names(forig))[di]
    rti <- grep("^rt", names(forig))[di]
    Qi <- grep("^Q", names(forig))[di]
    adducti <- grep("^adduct", names(forig))[di]
    cTable[lost, fields] <- forig[,c(idi, mzi, rti, Qi, adducti)]
    cTable[lost, "group"] <- as.integer(0)
    if("subgroup" %in% names(cTable))
        cTable[lost,"subgroup"] <- as.integer(0)
    samps <- unlist(getSamples(object)[datasets[[type]]])
    extras <- unlist(getExtra(object)[datasets[[type]]])
    samples_extras <- setdiff(c(samps, extras), combinerNames())
    if(any(duplicated(samples_extras))){
        warning("duplicate column names detected; information for these columns
                not copied as this could lead to erroneous behavior")
        dupnames <- names(which(table(samples_extras) > 1))
        samples_extras <- setdiff(samples_extras, dupnames)
    }
    samples_extras <- intersect(samples_extras, names(origTable))
    cTable[lost, samples_extras] <- origTable[samples_extras]

    return(cTable)
}



#' Update combinedTable & featdata
#'
#' @param cTable combinedTable data frame of the current metabCombiner object
#'
#' @param fdata featdata of the current metabCombiner object
#'
#' @param origdata xdata or ydata used to update the metabCombiner object
#'
#' @param datasets dataset identifiers in list format (x or y)
#'
#' @param type either one of 'x' (for xdata) or 'y' (for ydata)
#'
#' @noRd
updateCombiner <- function(object, cTable, fdata, origdata, datasets, type)
{
    ###object checks
    if(methods::is(origdata, "metabCombiner")){
        combinerCheck(isMetabCombiner(origdata), "metabCombiner")
        n <- length(datasets[[type]])
        if(length(datasets(origdata)) != n)
            stop("expected ", n, " datasets in input argument '",
                paste(type,"data", sep = ""),"'; found ",
                length(datasets(origdata)))
        forig <- featdata(origdata)
        origTable <- combinedTable(origdata)
    }

    else if(methods::is(origdata, "metabData")){
        combinerCheck(isMetabData(origdata), "metabData")
        n <- 1
        forig <- getData(origdata)[c("id","mz","rt","Q", "adduct")]
        origTable <-within(getData(origdata),rm("id","mz","rt","Q","adduct"))
    }

    if(type == "x")
        di <- seq(1,n)
    else if(type == "y")
        di <- seq(length(datasets$x)+1, length(datasets$x)+n)

    lostfeats <- find_lost_features(fdata, forig, type, di)
    fdata[seq(nrow(fdata)+1, nrow(fdata) + length(lostfeats)),] <- NA
    cTable[seq(nrow(cTable)+1, nrow(cTable) + length(lostfeats)),] <- NA
    fdata <- fill_featdata(fdata, forig, lostfeats, type, di)
    cTable <- fill_combinedTable(cTable, origTable, forig, lostfeats, datasets,
                                 type, object)
    return(list(fdata = fdata, cTable = cTable))
}



#' @title Update \code{metabCombiner} Objects
#'
#' @description This method updates the feature list (featdata) and aligned
#' table (\code{combinedTable}) within a \code{metabCombiner} object. Manual
#' changes to the (\code{combinedTable}) as well as unmatched X & Y dataset
#' features can be incorporated into the object and the corresponding results.
#' This function is typically paired with \code{link{reduceTable}} or other
#' forms of table reduction performed by the user.
#'
#' @param object \code{metabCombiner} object to be updated
#'
#' @param xdata \code{metabData} or \code{metabCombiner} object originally used
#' to construct the object argument
#'
#' @param ydata \code{metabData} or \code{metabCombiner} object originally used
#' to construct the object argument
#'
#' @param combinedTable merged table which may be altered by the user. This must
#' have the \code{combinedTable} format to be valid (see: ?isCombinedTable)
#'
#' @return \code{metabCombiner} object with updates to \code{combinedTable} to
#' include features that have been missed or changes by the user.
#'
#' @details There are two points where features can be removed from the
#' \code{combinedTable} report: during m/z grouping and during the table
#' reduction step. It is also possible for user-specified changes to the report
#' to remove certain features entirely. This function allows for the missed
#' features to be brought back into the table as non-matched entities. For xdata
#' features, the Y columns will be entirely missing values, and ydata features
#' will have missing X information. The feature data (featdata) will also be
#' updated for use in subsequent alignments, but only features present in the
#' representative dataset will be retained by default.
#'
#' @note Duplicated sample & extra column names cannot be copied from the
#' original data they feature in, therefore they are left as missing values.
#'
#' @examples
#' data(plasma30)
#' data(plasma20)
#' p30 <- metabData(plasma30, samples = "CHEAR")
#' p20 <- metabData(plasma20, samples = "Red")
#' p.comb <- metabCombiner(xdata = p30, ydata = p20, xid = "p30", yid = "p20",
#'                         binGap = 0.0075)
#'
#' ##extracting, modifying, and updating combinedTable
#' cTable <- combinedTable(p.comb)
#' cTable <- dplyr::filter(cTable, rty < 17.25)
#' p.comb <- updateTables(p.comb, combinedTable = cTable)
#'
#' p.comb <- selectAnchors(p.comb, tolmz = 0.003, tolQ = 0.3, windy = 0.02)
#' p.comb <- fit_gam(p.comb, k = 20, iterFilter = 1)
#' p.comb <- calcScores(p.comb, A = 90, B = 14, C = 0.5)
#' p.comb <- reduceTable(p.comb, delta = 0.2, maxRTerr = 0.5)
#'
#' ##updating to include features removed from xdata & ydata
#' p.comb <- updateTables(p.comb, xdata = p30, ydata = p20)
#'
#' #view results
#' cTable <- combinedTable(p.comb)
#' fdata <- featdata(p.comb)
#'
#' @export
updateTables <- function(object, xdata = NULL, ydata = NULL,
                         combinedTable = NULL)
{
    combinerCheck(isMetabCombiner(object), "metabCombiner")
    if(is.null(combinedTable)) cTable <- combinedTable(object)
    else cTable <- combinedTable
    fdata <- featdata(object)
    xy <- xy(object)
    dats <- datasets(object, list = TRUE)

    if(!is.null(combinedTable)){
        combinerCheck(isCombinedTable(cTable), "combinedTable")
        if(length(setdiff(cTable[["rowID"]], fdata[["rowID"]])) > 0)
            stop("at least one invalid rowID in combinedTable input")
        fdata <- fdata[match(cTable[["rowID"]], fdata[["rowID"]]),]
    }

    if(!is.null(xdata)){
        if(!methods::is(xdata, xy[["xtype"]]))
            stop("expected ", xy[["xtype"]], " for 'xdata' argument")
        updated_with_X <- updateCombiner(object,cTable,fdata,xdata, dats,"x")
        cTable <- updated_with_X[["cTable"]]
        fdata <- updated_with_X[["fdata"]]
    }


    if(!is.null(ydata)){
        if(!methods::is(ydata, xy[["ytype"]]))
            stop("expected ", xy[["ytype"]], " for 'ydata' argument")
        updated_with_Y <- updateCombiner(object,cTable,fdata, ydata,dats, "y")
        cTable <- updated_with_Y[["cTable"]]
        fdata <- updated_with_Y[["fdata"]]
    }

    object <- update_mc(object, combinedTable = cTable, featdata = fdata)
    return(object)
}


