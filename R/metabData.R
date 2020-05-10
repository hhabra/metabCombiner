#' @title Constructor for the metabData object.
#'
#' @description
#' This is a constructor for objects of type \code{metabData}.
#'
#' @param table  Path to file containing feature table or data.frame object
#'               containing features
#'
#' @param mz    Character name(s) or regular expression associated with data
#' column containing m/z values. The first column whose name contains this
#' expression will be selected for analysis.
#'
#' @param rt    Character name(s) or regular expression associated with data
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
#' keywords given, will search for longest stretch of numeric columns remaining.
#'
#' @param extra  Character names of columns containing additional feature
#' information, e.g.  non-analyzed sample values. All columns containing these
#' keywords are selected for analysis.
#'
#' @param Q   Character name(s) or regular expression associated with numeric
#' feature abundance quantiles. If NULL, abundance quantiles will be calculated
#' from sample intensities.
#'
#' @param rtmin     Numeric. Minimum retention time for analysis.
#'
#' @param rtmax    Numeric. Maximum retention time for analysis.
#'
#' @param misspc     Numeric. Threshold missingness percentage for analysis.
#'
#' @param zero      Logical. Whether to consider zero values as missing.
#'
#' @param measure Character. Central abundance measure, either "median" or "mean".
#'
#' @param duplicate  Numeric ordered pair (m/z, rt) feature tolerances. Pairs of
#' features within these tolerances are deemed duplicates and one of the pair is
#' removed (see: \code{\link{findDuplicates}})
#'
#' @details
#' A processed metabolomics feature data frame must contain columns for m/z, rt,
#' and numeric sample intensities. The user may also supply some optional fields
#' such as identity \code{id} and \code{adduct} label columns. Non-analyzed columns
#' can be included into the final output by specifying the names of these columns
#' in the \code{extra} argument. All required arguments are checked for validity
#' (e.g. no negative m/z or rt values, each column is used at most once,
#' column types are valid, etc...).
#'
#' Following construction is a pre-analysis filtering of rows that are either,
#' 1) Outside of the specified retention time range (\code{rtmin},\code{rtmax}),
#' 2) Missing in excess of \code{misspc} percent of analyzed samples, or
#' 3) deemed duplicates by small pairwise differences (mzTol,rtTol) as specified
#' by the \code{duplicate} argument. Remaining features are ranked by abundance
#' quantiles, Q, using a central \code{measure}, either "median" or "mean."
#' Alternatively, the abundance quantiles column can be specified in the
#' argument \code{Q}.
#'
#' @return An object of class metabData containing the specific information
#' specified by \code{mz,rt, samples, id, adduct, Q, and extra} arguments, and
#' adjusted by pre-processing steps.
#'
#' @examples
#' \dontrun{
#' library(metabCombiner)
#'
#' data(plasma30)
#' data(plasma20)
#'
#' #analyzing CHEAR plasma samples; RedCross samples non-analyzed "extra" columns
#' p30 <- metabData(plasma30, mz = "mz", rt = "rt", id = "identity", adduct = "adduct",
#'                  samples = "CHEAR", extra = "RedCross")
#'
#' #equivalent
#' p30 <- metabData(plasma30, id = "id", samples = "CHEAR", extra = "Red")
#'
#'
#' getSamples(p30)  #should print names of 5 CHEAR Sample column names
#' getExtra(p30)    #should print names of 5 Red Cross Sample column names
#'
#' #analyzing Red Cross samples with retention time limitations (0.5-17.5min)
#' p20 <- metabData(plasma20, samples = "Red", rtmin = 0.5, rtmax = 17.5)
#' data = getData(p20)
#' range(data$rt)
#'
#' #using regular expressions for field searches
#' p30.2 <- metabData(plasma30, id = "identity|id|ID", samples = ".[3-5]$")
#' getSamples(p30.2)    #should print all column names ending in .3, .4, .5
#' }
#' @export
metabData <- function(table, mz = "mz", rt = "rt", id = "id",
                      adduct = "adduct", samples = NULL, Q = NULL,
                      extra = NULL, rtmin = "min", rtmax = "max",
                      misspc = 50, measure = c("median", "mean"),
                      zero = FALSE, duplicate = c(0.0025, 0.05))
{
    if(missing(table))
        stop("required argument 'table' is missing with no default")

    if(!is.character(mz))
        stop("non-character argument for variable 'mz'")

    if(!is.character(rt))
        stop("non-character argument for variable 'rt'")

    if(!is.character(samples))
        stop("non-character argument for variable 'samples'")

    if(!is.character(id))
        id = NULL

    if(!is.character(adduct))
        adduct = NULL

    if(!is.character(extra))
        extra = NULL

    if (misspc >= 100 | misspc < 0 | !is.numeric(misspc))
        stop("Parameter 'misspc' must be a numeric value from [0,100)")

    if(typeof(table) == "character")
        table <- readData(table)
    else if(dplyr::is.tbl(table))     #handling tbl error
        table <- as.data.frame(table)
    else if(!is.data.frame(table))
        stop("argument 'table' must be a data.frame or path to data file")

    if(!is.logical(zero)){
        warning("non-logical value for argument 'zero'. Setting to default (FALSE)")
        zero = FALSE
    }

    measure = match.arg(measure)

    newData <- new("metabData")

    newData = detectFields(Data = newData, table = table,
                           mz = mz, rt = rt, id, adduct = adduct,
                           samples = samples, extra = extra, Q = Q)

    newData = adjustData(Data = newData, misspc = misspc, measure = measure,
                         rtmin = rtmin, rtmax = rtmax, zero = zero,
                         duplicate = duplicate)

    return(newData)
}






