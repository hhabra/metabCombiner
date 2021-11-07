############  Generics put in this file in alphabetical order #############

#' @title Retrieve Adduct Annotations
#'
#' @description  This retrieves user-assigned adduct annotations from one or
#' all constituent datasets of a \code{metabCombiner} object
#'
#' @param object \code{metabCombiner} object
#'
#' @param data dataset identifier to extract information from; if NULL, extracts
#' information from all datasets
#'
#' @return data frame of adduct annotations
#'
#' @examples
#' data(plasma30)
#' data(plasma20)
#'
#' p30 <- metabData(head(plasma30,500), samples = "CHEAR")
#' p20 <- metabData(head(plasma20,500), samples = "Red")
#'
#' p.comb <- metabCombiner(p30, p20, xid = "p30", yid = "p20")
#'
#' ##retrieve all adduct data
#' adducts <- adductdata(p.comb, data = NULL)
#'
#' ##retrieve adduct data from p30
#' adducts <- adductdata(p.comb, data = "p30")
#'
#' @export
setGeneric("adductdata", function(object, data = NULL)
            standardGeneric("adductdata"))


#' @title Obtain Feature Alignment Report
#'
#' @description Obtain constructed table reporting every possible metabolomics
#' feature pair alignment.
#'
#' @param object \code{metabCombiner} object.
#'
#' @return Feature Pair Alignment report data frame. The columns of the report
#' are as follows:
#'
#' \item{idx}{Identities of features from dataset X}
#' \item{idy}{Identities of features from dataset Y}
#' \item{mzx}{m/z values of features from dataset X}
#' \item{mzy}{m/z values of features from dataset Y}
#' \item{rtx}{retention time values of features from dataset X}
#' \item{rty}{retention time values of features from dataset Y}
#' \item{rtProj}{model-projected (X->Y) retention times values}
#' \item{Qx}{abundance quantile values of features from dataset X}
#' \item{Qy}{abundance quantile values of features from dataset Y}
#' \item{group}{m/z feature group of feature pairing}
#' \item{score}{computed similarity scores of feature pairing}
#' \item{rankX}{ranking of pairing score for X dataset features}
#' \item{rankY}{ranking of pairing score for Y dataset features}
#' \item{adductX}{adduct label of features from dataset X}
#' \item{adductY}{adduct label of features from dataset Y}
#' \item{...}{Sample and extra columns from both datasets X & Y}
#'
#' @examples
#' data(plasma30)
#' data(plasma20)
#'
#' p30 <- metabData(head(plasma30,500), samples = "CHEAR")
#' p20 <- metabData(head(plasma20,500), samples = "Red")
#'
#' p.comb <- metabCombiner(p30, p20)
#' p.comb.table <- combinedTable(p.comb)
#'
#' @export
setGeneric("combinedTable", function(object) standardGeneric("combinedTable"))


#' @title Obtain Dataset IDs
#'
#' @description Each dataset in a \code{metabCombiner} object is represented by
#' a character identifier. The datasets slot contains all these ids in a single
#' vector, which can be obtained in sequential order with this accessor method
#'
#' @param object metabCombiner object
#'
#' @param list logical, option to return in list format (TRUE) vs character
#' vector format (FALSE)
#'
#' @return character vector of dataset identifiers
#'
#' @examples
#' #' @examples
#' data(plasma30)
#' data(plasma20)
#'
#' p30 <- metabData(head(plasma30,500), samples = "CHEAR")
#' p20 <- metabData(head(plasma20,500), samples = "Red")
#'
#' p.comb <- metabCombiner(p30, p20, xid = "p30", yid = "p20")
#'
#' ##datasets extraction: expect "p30", "p20"
#' sets <- datasets(p.comb, list = FALSE)
#'
#' @export
setGeneric("datasets", function(object, list = FALSE)
    standardGeneric("datasets"))


#' @title Obtain Feature Metadata
#'
#' @description \code{metabCombiner} objects organize metabolomics feature
#' information in the "featdata" slot. This method retrieves all metadata
#' or that of one dataset. The rows should identically correspond to the same
#' rows from the combinedTable data frame.
#'
#' @param object a \code{metabCombiner} object
#'
#' @param data character dataset identifier
#'
#' @return data frame of feature metadata from one or all datasets
#'
#' @examples
#' data(plasma30)
#' data(plasma20)
#'
#' p30 <- metabData(head(plasma30,500), samples = "CHEAR")
#' p20 <- metabData(head(plasma20,500), samples = "Red")
#'
#' p.comb <- metabCombiner(p30, p20, xid = "p30", yid = "p20")
#'
#' #full metadata extraction
#' fdata <- featdata(p.comb, data = NULL)
#'
#' #single dataset feature information extraction
#' fdata <- featdata(p.comb, data = "p20")
#'
#'@export
setGeneric("featdata",function(object, data = NULL) standardGeneric("featdata"))


#' @title Get Ordered Retention Time Pairs
#'
#' @description
#' This returns the data frame of feature alignments used to anchor a retention
#' time projection model, constructed by \code{\link{selectAnchors}}.
#'
#' @param object metabCombiner object
#'
#' @return Data frame of anchor features
#'
#' @seealso
#' \code{\link{selectAnchors}}
#'
#' @examples
#' data(plasma30)
#' data(plasma20)
#'
#' p30 <- metabData(plasma30, samples = "CHEAR")
#' p20 <- metabData(plasma20, samples = "Red")
#'
#' p.comb <- metabCombiner(p30, p20)
#' p.comb <- selectAnchors(p.comb, windx = 0.05, windy = 0.03)
#'
#' anchors <- getAnchors(p.comb)
#'
#' @export
setGeneric("getAnchors", function(object) standardGeneric("getAnchors"))


#' @title Obtain Last-Used Score Coefficients
#'
#' @description
#' Provides the last used weight arguments from \code{calcScores()} function.
#' Returns empty list if \code{calcScores()} has not yet been called.
#'
#' @param object  metabCombiner object
#'
#' @return A list of the last used weight parameters:
#' \item{A}{Specific weight penalizing feature m/z differences}
#' \item{B}{Specific weight penalizing retention time projection error}
#' \item{C}{Specific weight penalizing differences in abundance quantiles}
#'
#' @examples
#' data(plasma30)
#' data(plasma20)
#'
#' p30 <- metabData(plasma30, samples = "CHEAR")
#' p20 <- metabData(plasma20, samples = "Red")
#'
#' p.comb <- metabCombiner(p30, p20)
#' p.comb <- selectAnchors(p.comb, windx = 0.05, windy = 0.04, tolrtq = 0.15)
#' p.comb <- fit_gam(p.comb, k = 20, iterFilter = 1, family = "gaussian")
#' p.comb <- calcScores(p.comb, A = 90, B = 14, C = 0.5)
#'
#' getCoefficients(p.comb)
#'
#' @export
setGeneric("getCoefficients", function(object)
            standardGeneric("getCoefficients"))

#' Get Processed Dataset
#'
#' @description The \code{\link{metabData}} constructor creates a formatted
#' dataset from the input, which may be accessed using this method.
#'
#' @param object metabData object
#'
#' @return Single Metabolomics Data Frame
#'
#' @examples
#' data(plasma30)
#'
#' p30 <- metabData(plasma30, samples = "CHEAR")
#' data <- getData(p30)
#'
#' @export
setGeneric("getData", function(object) standardGeneric("getData"))


#' @title Get Fitted RT Model
#'
#' @description
#' Returns the last fitted RT projection model from a metabCombiner object of
#' type "gam" or "loess".
#'
#' @param object  metabCombiner object
#'
#' @param fit   Choice of model, "gam" or "loess"
#'
#' @return nonlinear retention time fit object
#'
#' @seealso
#' \code{\link{fit_gam}}, \code{\link{fit_loess}}
#'
#' @examples
#' data(plasma30)
#' data(plasma20)
#' p30 <- metabData(plasma30, samples = "CHEAR")
#' p20 <- metabData(plasma20, samples = "Red", rtmax = 17.25)
#' p.comb <- metabCombiner(xdata = p30, ydata = p20, binGap = 0.005)
#' p.comb <- selectAnchors(p.comb, tolrtq = 0.15, tolQ = 0.2, windy = 0.02)
#' p.comb <- fit_gam(p.comb, iterFilter = 1, k = 20, family = "gaussian")
#' p.comb <- fit_loess(p.comb, iterFilter = 1, spans = 0.2)
#' model.gam <- getModel(p.comb, fit = "gam")
#' model.loess <- getModel(p.comb, fit = "loess")
#'
#' @export
setGeneric("getModel", function(object, fit = c("gam", "loess"))
            standardGeneric("getModel"))


#' Get Extra Data Column Names
#'
#' @param object  \code{metabCombiner} or \code{metabData} object
#'
#' @param data dataset identifier for \code{metabCombiner} objects
#'
#' @return character vector of extra column names
#'
#' @examples
#' data(plasma30)
#' p30 <- metabData(plasma30, samples = "CHEAR", extra = "Red")
#' getExtra(p30)
#'
#' @export
setGeneric("getExtra", function(object, data = NULL)
            standardGeneric("getExtra"))


#' @title Get Sample Names From metabCombiner or metabData Object
#'
#' @description Returns the sample names from one of the two datasets used in
#' metabCombiner analysis, denoted as 'x' or 'y.'
#'
#' @param object  \code{metabCombiner} or \code{metabData} object
#'
#' @param data dataset identifier for \code{metabCombiner} objects
#'
#' @return character vector of sample names. For \code{metabCombiner} objects
#'     these may come from the 'x' dataset (if \code{data} = "x") or the 'y'
#'     dataset (if \code{data} = "y").
#'
#' @examples
#' data(plasma30)
#' data(plasma20)
#'
#' p30 <- metabData(plasma30, samples = "CHEAR")
#' p20 <- metabData(plasma20, samples = "Red", rtmax = 17.25)
#'
#' p.comb <- metabCombiner(xdata = p30, ydata = p20)
#'
#' getSamples(p30)
#' getSamples(p.comb, data = "x")  #equivalent to previous
#' getSamples(p20)
#' getSamples(p.comb, data = "y")  #equivalent to previous
#'
#' @export
setGeneric("getSamples", function(object, data = NULL)
            standardGeneric("getSamples"))


#' Get Object Statistics
#'
#' @description Prints out a list of object-specific statistics for both
#' \code{metabCombiner} and \code{metabData} objects
#'
#' @param object metabCombiner or metabData object
#'
#' @return list of object-specific statistics
#'
#' @examples
#' data(plasma30)
#' data(plasma20)
#' p30 <- metabData(plasma30, samples = "CHEAR")
#' p20 <- metabData(plasma20, samples = "Red", rtmax = 17.25)
#'
#' getStats(p30) #metabData stats
#'
#' p.comb <- metabCombiner(xdata = p30, ydata = p20, binGap = 0.005)
#' p.comb <- selectAnchors(p.comb, tolmz = 0.003, tolQ = 0.3, windy = 0.02)
#' p.comb <- fit_gam(p.comb, iterFilter = 1, k = 20)
#'
#' getStats(p.comb) #metabCombiner stats
#'
#' @export
setGeneric("getStats", function(object) standardGeneric("getStats"))


#' @title Retrieve Feature Identities
#'
#' @description  This retrieves user-assigned feature identities from one or
#' all constituent datasets of a \code{metabCombiner} object
#'
#' @param object \code{metabCombiner} object
#'
#' @param data dataset identifier to extract information from; if NULL, extracts
#' information from all datasets
#'
#' @return data frame of feature identities
#'
#' @examples
#' data(plasma30)
#' data(plasma20)
#'
#' p30 <- metabData(head(plasma30,500), samples = "CHEAR")
#' p20 <- metabData(head(plasma20,500), samples = "Red")
#' p.comb <- metabCombiner(p30, p20, xid = "p30", yid = "p20")
#'
#' ##retrieve all ids
#' ids <- iddata(p.comb, data = NULL)
#'
#' ##retrieve ids from p30
#' ids <- iddata(p.comb, data = "p30")
#'
#' @export
setGeneric("iddata", function(object, data = NULL)
    standardGeneric("iddata"))


#' @title Retrieve m/z Values
#'
#' @description  This retrieves feature m/z values from one or all constituent
#' datasets of a \code{metabCombiner} object. Alternatively, the average
#' m/z value can be retrieved.
#'
#' @param object \code{metabCombiner} object
#'
#' @param data dataset identifier to extract information from; if NULL,
#' extracts data frame information from all datasets
#'
#' @param value Either "obs" (observed - default option) or "mean" value
#'
#' @return data frame of m/z values (if NULL) or single vector of m/z values
#'
#' data(plasma30)
#' data(plasma20)
#'
#' p30 <- metabData(head(plasma30,500), samples = "CHEAR")
#' p20 <- metabData(head(plasma20,500), samples = "Red")
#' p.comb <- metabCombiner(p30, p20, xid = "p30", yid = "p20")
#'
#' ##retrieve all m/z
#' mz <- mzdata(p.comb, data = NULL)
#'
#' ##retrieve m/z from p30
#' mz <- mzdata(p.comb, data = "p30")
#'
#' ##retrieve mean m/z
#' mz <- mzdata(p.comb, value = "mean")
#'
#' @export
setGeneric("mzdata", function(object, data = NULL, value = c("obs", "mean"))
    standardGeneric("mzdata"))


#' @title Get Nonmatched Features
#'
#' @description
#' Features that lack a any counterparts in the complementary dataset may be
#' obtained from this method. If data is set to "x" or "y", will retrieve data
#' from the current X or Y dataset, respectively. If data is set to NULL, will
#' retrieve the list of nonmatched features.
#'
#' @param object  metabCombiner object
#'
#' @param data dataset identifier for \code{metabCombiner} objects; if NULL,
#' returns full list of non-matched features
#'
#' @return Data frame of non-matched features corresponding to data argument
#'
#' @examples
#' data(plasma30)
#' data(plasma20)
#'
#' p30 <- metabData(head(plasma30,500), samples = "CHEAR")
#' p20 <- metabData(head(plasma20,500), samples = "Red", rtmax = 17.25)
#' p.comb <- metabCombiner(xdata = p30, ydata = p20, binGap = 0.005)
#'
#' nnmx <- nonmatched(p.comb, data = "x")
#' nnmy <- nonmatched(p.comb, data = "y")
#'
#' @export
setGeneric("nonmatched", function(object, data = "x")
            standardGeneric("nonmatched"))


#' @title Retrieve Relative Abundance Values
#'
#' @description  This retrieves feature Q values from one or all constituent
#' dataset features of a \code{metabCombiner} object. Alternatively, the average
#' Q value can be retrieved.
#'
#' @param object \code{metabCombiner} object
#'
#' @param data dataset identifier to extract information from; if NULL,
#' extracts information from all datasets
#'
#' @param value Either "obs" (observed - default option) or "mean" average value
#'
#' @return data frame or vector of relative ranked abundance (Q) values
#'
#' data(plasma30)
#' data(plasma20)
#'
#' p30 <- metabData(head(plasma30,500), samples = "CHEAR")
#' p20 <- metabData(head(plasma20,500), samples = "Red")
#' p.comb <- metabCombiner(p30, p20, xid = "p30", yid = "p20")
#'
#' ##retrieve all Q
#' Q <- Qdata(p.comb, data = NULL)
#'
#' ##retrieve Q from p30
#' Q <- Qdata(p.comb, data = "p30")
#'
#' ##retrieve mean Q
#' Q <- Qdata(p.comb, value = "mean")
#'
#' @export
setGeneric("Qdata", function(object, data = NULL, value = c("obs", "mean"))
    standardGeneric("Qdata"))


#' @title Retrieve Retention Time Values
#'
#' @description  This retrieves feature RT values from one or all constituent
#' dataset features of a \code{metabCombiner} object. Alternatively, the average
#' RT value can be retrieved.
#'
#' @param object \code{metabCombiner} object
#'
#' @param data dataset identifier to extract information from; if NULL,
#' extracts information from all datasets
#'
#' @param value Either"obs" (observed - default option) or "mean"
#'
#' @return data frame or vector of retention time values
#'
#' @examples
#' data(plasma30)
#' data(plasma20)
#'
#' p30 <- metabData(head(plasma30,500), samples = "CHEAR")
#' p20 <- metabData(head(plasma20,500), samples = "Red")
#' p.comb <- metabCombiner(p30, p20, xid = "p30", yid = "p20")
#'
#' ##retrieve all RTs
#' rt <- rtdata(p.comb, data = NULL)
#'
#' ##retrieve RTs from p30
#' rt <- rtdata(p.comb, data = "p30")
#'
#' ##retrieve mean RT
#' rt <- rtdata(p.comb, value = "mean")
#'
#' @export
setGeneric("rtdata", function(object, data = NULL, value = c("obs", "mean"))
    standardGeneric("rtdata"))


##updater methods
setGeneric("update_mc",function(object, combinedTable, featdata, anchors, model,
                          fit, samples, extra, xy, datasets, nonmatched, stats,
                          values, coefficients) standardGeneric("update_mc"))


setGeneric("update_md", function(object, data, samples, extra, stats)
            standardGeneric("update_md"))



#' @title Obtain x & y Data Identifiers
#'
#' @description \code{metabCombiner} alignment is performed in a pairwise manner
#' between two datasets generically termed "x" & "y". These methods print the
#' identifier(s) associated with datasets X and Y, contained within the xy slot
#' of a constructed \code{metabCombiner} object.
#'
#' @param object \code{metabCombiner} object
#'
#' @return character X or Y dataset identifiers
#'
#' data(plasma30)
#' data(plasma20)
#'
#' p30 <- metabData(head(plasma30,500), samples = "CHEAR")
#' p20 <- metabData(head(plasma20,500), samples = "Red")
#' p.comb <- metabCombiner(p30, p20, xid = "p30", yid = "p20")
#'
#' #expected: "p30"
#' x(p.comb)
#'
#' #expected: "p20"
#' y(p.comb)
#'
#' #list of x & y data descriptors
#' xy(p.comb)
#'
#' @export
setGeneric("x", function(object) standardGeneric("x"))

#' @rdname x
#'
#' @export
setGeneric("xy", function(object) standardGeneric("xy"))

#' @rdname x
#'
#' @export
setGeneric("y", function(object) standardGeneric("y"))
