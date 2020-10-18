#' Combiner: Method for Combining LC-MS Metabolomics Feature Measurements
#'
#' Takes two untargeted metabolomics feature lists (consisting of m/z,
#' rt, and sample intensity measurements, plus optional identifiers & adduct
#' labels) and outputs a merged feature list consisting of potential
#' compound matches, ranked by a similarity score for groups of features. Inputs
#' are assumed to be derived from biologically similar samples analyzed with
#' a similar analytical method.
#'
#' @docType package
#' @name metabCombiner
#'
#' @useDynLib metabCombiner
#' @import dplyr
#' @importFrom graphics plot lines
#' @importFrom stats loess loess.control predict setNames
#' @importFrom S4Vectors isSorted
#' @importFrom methods new
#' @importFrom utils read.csv read.delim write.table
#' @importFrom caret createFolds
#' @importFrom rlang .data
#' @importFrom matrixStats rowCounts rowMedians rowMeans2
#'
NULL

