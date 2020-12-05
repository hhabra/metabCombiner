############  Generics put in this file in alphabetical order #############

#' @title Obtain metabCombiner Feature Alignment Report
#'
#' @description Obtain constructed table reporting every possible metabolomics
#'     feature alignment.
#'
#' @param object metabCombiner object.
#'
#' @return Feature Pair Alignment report data.frame. The columns of the report
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
#' p30 <- metabData(plasma30, samples = "CHEAR")
#' p20 <- metabData(plasma20, samples = "Red")
#'
#' p.comb <- metabCombiner(p30, p20)
#' p.comb.table <- combinedTable(p.comb)
#'
#' @export
setGeneric("combinedTable", function(object) standardGeneric("combinedTable"))


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
#' p.comb <- selectAnchors(p.comb, windx = 0.05, windy = 0.03)
#' p.comb <- fit_gam(p.comb, k = 20, iterFilter = 1)
#' p.comb <- calcScores(p.comb, A = 90, B = 14, C = 0.5)
#'
#' getCoefficients(p.comb)
#'
#' @export
setGeneric("getCoefficients", function(object)
            standardGeneric("getCoefficients"))

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
#' p.comb <- selectAnchors(p.comb, tolmz = 0.003, tolQ = 0.3, windy = 0.02)
#' p.comb <- fit_gam(p.comb, iterFilter = 1, k = 20)
#' p.comb <- fit_loess(p.comb, iterFilter = 1, spans = 0.2)
#' model.gam <- getModel(p.comb, fit = "gam")
#' model.loess <- getModel(p.comb, fit = "loess")
#'
#' @export
setGeneric("getModel", function(object, fit = c("gam", "loess"))
            standardGeneric("getModel"))

#' Get Processed Dataset
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


#' Get Extra Data Column Names
#'
#' @param object metabCombiner or metabData object
#'
#' @param data Choice of input dataset, e.g. "x" or "y"
#'
#' @return character vector of extra column names
#'
#' @examples
#' data(plasma30)
#' p30 <- metabData(plasma30, samples = "CHEAR", extra = "Red")
#' getExtra(p30)
#'
#' @export
setGeneric("getExtra", function(object, data = "x")
            standardGeneric("getExtra"))


#' @title Get Sample Names From metabCombiner or metabData Object
#'
#' @description Returns the sample names from one of the two datasets used in
#' metabCombiner analysis, denoted as 'x' or 'y.'
#'
#' @param object  metabCombiner or metabData object
#'
#' @param data   Character. One of either 'x' or 'y'.
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
setGeneric("getSamples", function(object, data = "x")
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
#'
#' @export
setGeneric("getStats", function(object) standardGeneric("getStats"))

#' @title Get Nonmatched Features
#'
#' @description
#' Features that lack a any counterparts in the complementary dataset may be
#' obtained from this method.
#'
#' @param object  metabCombiner object
#'
#' @param data  Either one of 'x' or 'y', specifying which dataset's nonmatched
#'      features to return.
#'
#' @return If data is "x", returns non-matched X features ; if "y", returns
#'         non-matched Y features
#'
#' @examples
#' data(plasma30)
#' data(plasma20)
#'
#' p30 <- metabData(plasma30, samples = "CHEAR")
#' p20 <- metabData(plasma20, samples = "Red", rtmax = 17.25)
#' p.comb <- metabCombiner(xdata = p30, ydata = p20, binGap = 0.005)
#'
#' nnmx <- nonmatched(p.comb, data = "x")
#' nnmy <- nonmatched(p.comb, data = "y")
#'
#' @export
setGeneric("nonmatched", function(object, data = "x")
            standardGeneric("nonmatched"))

##updater methods
setGeneric("update_mc",function(object, data, combinedTable, nonmatched,
            samples,extra,anchors, fit, model, coefficients, stats, values)
            standardGeneric("update_mc"))


setGeneric("update_md", function(object, data, samples, extra, stats)
            standardGeneric("update_md"))


