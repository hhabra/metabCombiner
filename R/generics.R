############  Generics put in this file in alphabetical order #############

#' @title Obtain metabCombiner Feature Alignment Report
#'
#' @description Obtain constructed table reporting every possible metabolomics
#' feature alignment.
#'
#' @param object metabCombiner object.
#'
#' @return  Feature Pair Alignment report data.frame. The columns of the report
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
#' \item{B}{Specific weight penalizing relative error of retention time projection}
#' \item{C}{Specific weight penalizing differences in abundance quantiles}
#'
#' @export
setGeneric("getCoefficients", function(object) standardGeneric("getCoefficients"))

#' @title Get Fitted RT Model
#'
#' @description
#' Returns the last fitted RT projected model from a metabCombiner object of
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
#' @export
setGeneric("getModel", function(object, fit = c("gam", "loess")) standardGeneric("getModel"))


#' Get Processed Dataset
#'
#' @param object metabData object
#'
#' @return Single Metabolomics Data Frame
#'
#' @export
setGeneric("getData", function(object) standardGeneric("getData"))


#' Get Extra Data Column Names
#'
#' @param object metabCombiner or metabData object
#'
#' @param data Choice of input dataset, 'x' or 'y'
#'
#' @return character vector of extra column names
#'
#' @export
setGeneric("getExtra", function(object, data = c("x", "y")) standardGeneric("getExtra"))


#' @title Get Sample Names From metabCombiner Object
#'
#' @description \code{metabCombiner} objects consist of two formatted metabolomics
#' feature tables. This method returns the sample names from one of the two datasets.
#'
#' @param object  metabCombiner or metabData object
#'
#' @param data   Character. One of either 'x' or 'y'.
#'
#' @return If data is "x", returns sample names for dataset X; if "y", returns
#' sample names from dataset Y.
#'
#' @export
setGeneric("getSamples", function(object, data = c("x", "y")) standardGeneric("getSamples"))


#' Get Object Statistics
#'
#' @description Prints out a list of object-specific statistics for both
#' \code{metabCombiner} and \code{metabData} objects
#'
#' @param object metabCombiner or metabData object
#'
#' @return list of object-specific statistics
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
#' @param data    Either one of 'x' or 'y', specifying which dataset's nonmatched
#'                features to return.
#'
#' @return If data is "x", returns non-matched X features ; if "y", returns
#'         non-matched Y features
#'
#' @export
setGeneric("nonmatched", function(object, data = c("x", "y")) standardGeneric("nonmatched"))
