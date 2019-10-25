##class descriptions go in here
#

############################# Classes #####################

#' @title 'metabData' Single Metabolomics Dataset Class
#'
#' @description
#' This class is designed to pre-process, format, and acilitating the subsequent
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
#' This is the main object for the \code{Combiner} package workflow. This
#' object holds a pair of analyzed metabolomics datasets and a combined
#' feature table, along with a retention time projection model.
#'
#' @slot xdata  One of two formatted metabData objects to be combined.
#'
#' @slot ydata  One of two formatted metabData objects to be combined.
#'
#' @slot combinerTable  A data.frame combining tables measurements and displaying
#' all feature alignments
#'
#' @slot binGap  The numeric m/z gap used to compute feature groups
#'
#' @slot anchors  A data.frame consisting of feature alignments to be used for
#' fitting a nonlinear retention time projection model
#'
#' @slot model  A list containing the last fitted nonlinear model(s).
#'
#' @slot coefficients A list of last used A,B,C similarity score weight arguments
#'
#' @slot stats A list of useful metabCombiner object statistics
#'
#' @export
setClass("metabCombiner", slots = c(xdata = "metabData",
                                    ydata = "metabData",
                                    combinerTable = "data.frame",
                                    binGap = "numeric",
                                    anchors = "data.frame",
                                    model = "list",
                                    coefficients = "list",
                                    stats = "list"
                          ),
                          prototype = prototype(
                                   xdata = new("metabData"),
                                   ydata = new("metabData"),
                                   combinerTable = data.frame(),
                                   binGap = numeric(),
                                   anchors = data.frame(),
                                   model = list(),
                                   coefficients = list(),
                                   stats = list()
                          )
)




