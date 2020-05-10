##class descriptions go in here

################################## Classes ################################

#' @title 'metabData' Single Metabolomics Dataset Class
#'
#' @description
#' This class is designed to pre-process, format, and facilitate the subsequent
#' combining of input metabolomics feature tables. It stores the information of
#' a single metabolomics dataset, including the formatted feature table, sample
#' names, and feature statistics.
#'
#' @slot data A formatted metabolomics data frame.
#'
#' @slot samples A character vector of analyzed samples.
#'
#' @slot extra A character vector of names pertaining to non-analyzed columns
#'
#' @slot stats A list of dataset statistics
#'
#' @export
setClass("metabData", slots = c(data = "data.frame",
                                samples = "character",
                                extra = "character",
                                stats = "list"),
                      prototype = prototype(
                                 data = data.frame(),
                                 samples = character(),
                                 extra = character(),
                                 stats = list()
                      )
)


#' @title 'metabCombiner' Combined Metabolomics Dataset Class
#'
#' @description
#' This is the main object for the \code{metabCombiner} package workflow. This
#' object holds a combined feature table, along with a retention time
#' projection model, the ordered pair anchors used to generate this model, and
#' key object statistics.
#'
#' @slot combinedTable  A data.frame combining tables measurements and displaying
#' all feature alignments
#'
#' @slot nongrouped A list of data frames containing features from input tables
#' that could not be paired in \code{combinedTable}
#'
#' @slot anchors  A data.frame consisting of feature alignments to be used for
#' fitting a nonlinear retention time projection model
#'
#' @slot model  A list containing the last fitted nonlinear model(s).
#'
#' @slot coefficients A list of last used A,B,C similarity weight values
#'
#' @slot samples A list of sample name character vectors from input datasets
#'
#' @slot extra A list of extra column chacter name vectors from input datasets
#'
#' @slot stats A set of useful metabCombiner statistics
#'
#' @export
setClass("metabCombiner", slots = c(combinedTable = "data.frame",
                                    nongrouped = "list",
                                    anchors = "data.frame",
                                    model = "list",
                                    coefficients = "list",
                                    samples = "list",
                                    extra = "list",
                                    stats = "list"
                          ),
                          prototype = prototype(
                                   combinedTable = data.frame(),
                                   nongrouped = list(),
                                   anchors = data.frame(),
                                   model = list(),
                                   coefficients = list(),
                                   samples = list(),
                                   extra = list(),
                                   stats = list()
                          )
)
