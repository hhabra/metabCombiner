##class descriptions go in here

################################## Classes ################################

#' @title 'metabData' Single Metabolomics Dataset Class
#'
#' @description
#' This class is designed to process and format input metabolomics feature
#' tables. It stores the information from individual metabolomics datasets,
#' including the formatted feature table, sample names, and feature statistics.
#'
#' @slot data formatted metabolomics data frame.
#'
#' @slot samples character vector of analyzed sample names
#'
#' @slot extra character vector of non-analyzed columns names
#'
#' @slot stats A list of dataset statistics
#'
#' @slot filtered A list of filtered dataset features
#'
#' @export
setClass("metabData", slots = c(data = "data.frame",
                                samples = "character",
                                extra = "character",
                                stats = "list",
                                filtered = "list"),
                    prototype = prototype(
                                data = data.frame(),
                                samples = character(),
                                extra = character(),
                                stats = list(),
                                filtered = list()
                    )
)


#' @title 'metabCombiner' Combined Metabolomics Dataset Class
#'
#' @description
#' This is the main object for the \code{metabCombiner} package workflow. This
#' object holds a combined feature table, along with a retention time warping
#' model, the ordered pair anchors used to generate this model, important
#' information organized by dataset, and key object statistics.
#'
#' @slot combinedTable data frame displaying all feature pair alignments,
#' combining measurements of all possible shared compounds
#'
#' @slot featData data frame of feature metadata (id, m/z, rt, Q, adduct)
#'
#' @slot anchors data frame of feature pairs used for RT warping model
#'
#' @slot model  list containing the last fitted nonlinear model(s)
#'
#' @slot datasets list of constituent datasets from xdata & ydata inputs
#'
#' @slot xy current X & Y datasets
#'
#' @slot nonmatched list of data frames consisting of nonmatched features
#'
#' @slot coefficients list of last used A,B,C similarity weight values
#'
#' @slot samples list of sample name vectors from input datasets
#'
#' @slot extra list of extra column name vectors from input datasets
#'
#' @slot stats set of useful metabCombiner statistics
#'
#' @export
setClass("metabCombiner", slots = c(combinedTable = "data.frame",
                                    featData = "data.frame",
                                    anchors = "data.frame",
                                    model = "list",
                                    datasets = "list",
                                    xy = "list",
                                    samples = "list",
                                    extra = "list",
                                    nonmatched = "list",
                                    coefficients = "list",
                                    stats = "list"
                            ),
                            prototype = prototype(
                                    combinedTable = data.frame(),
                                    featData = data.frame(),
                                    anchors = data.frame(),
                                    model = list(),
                                    datasets = list(),
                                    xy = list(),
                                    samples = list(),
                                    extra = list(),
                                    nonmatched = list(),
                                    coefficients = list(),
                                    stats = list()
                            )
)
