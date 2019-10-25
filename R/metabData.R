#' @title Constructor for the metabData object.
#'
#' @description
#' This is a constructor for objects of type \code{metabData}. A processed
#'   metabolomics feature data frame must contain columns for m/z, rt, and numeric
#'   sample intensities. The user may also supply some optional fields such as
#'   \code{id} and \code{adduct} label columns. Non-analyzed columns can be
#'   included into the final output by specifying the names of these columns in the
#'   "extra" argument. All required arguments are checked for correctness and
#'   validity (e.g. no negative m/z or rt values, no overlapping columns, etc...).
#'
#' Following construction is a pre-analysis filtering of rows that are either,
#'   1) Outside of the specified retention time range (\code{rtmin},\code{rtmax}),
#'   2) Missing in excess of \code{misspc} percent of analyzed samples, or 3)
#'   deemed duplicates by small pairwise differences (mzTol,rtTol) as specified
#'   by the \code{duplicate} argument. Remaining features are ranked by abundance
#'   quantiles, Q, using a central \code{measure}, either "median" or "mean."
#'   Alternatively, the abundance quantiles column can be specified in the
#'   argument \code{Q}.
#'
#' \code{metabData()} currently supports a limited set of options for missing value
#'   imputation, specified by the \code{impute} & \code{imputeVal} options.
#'
#' @param table  Path to file containing feature table or data.frame object
#'               containing features
#'
#' @param mz    Character. Name for column containing m/z values. If "mz",
#'              will search for {mz, m/z, m.z, mass}
#'
#' @param rt    Character. Name for column containing retention time values.
#'              If "rt", will search for {rt, retention time, r.t, RT}.
#'
#' @param id         Character. Name for column containing compound identifiers.
#'                  If "id", will search for {id, compound, feature}
#'
#' @param adduct     Character. Name for column containing adduct labels. If
#'                  "adduct", will search for {adduct, adducts, annotation}.
#'
#' @param samples  Character. Names of columns containing sample values. If
#'                "detect", finds longest stretch of consecutive numeric columns.
#'
#' @param Q     Character. Name of user supplied quantile column.
#'
#' @param extra      Character. Names of (additional) user-supplied columns.
#'
#' @param rtmin      Numeric. Minimum retention time for analysis.
#'
#' @param rtmax      Numeric. Maximum retention time for analysis.
#'
#' @param misspc     Numeric. Threshold missingness percentage for analysis.
#'
#' @param zero      Logical. Whether to consider zero values as missing.
#'
#' @param measure Character. Central abundance measure, either "median" or "mean".
#'
#' @param impute     Logical. Whether to impute value for missing sample measures.
#'
#' @param imputeVal  Imputed value. One of "median", "mean", or numeric value.
#'
#' @param duplicate  Numeric ordered pair (m/z, rt) duplicate feature tolerances.
#'
#' @return An object of class metabData containing the specific information
#' specified by \code{mz,rt, samples, id, adduct, Q, and extra} arguments, and adjusted
#' by pre-processing steps.
#'
#' @export
metabData <- function(table, mz = "mz", rt = "rt", id = "id",
                      adduct = "adduct", samples = NULL, Q = NULL,
                      extra = NULL, rtmin = "min", rtmax = "max",
                      misspc = 50, measure = c("median", "mean"),
                      zero = FALSE, impute = FALSE, imputeVal = 100,
                      duplicate = c(0.0025, 0.05))
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

    if(!is.logical(impute)){
        warning("Parameter 'impute' is non-logical. Setting to default (FALSE)")
        impute = FALSE
    }

    measure = match.arg(measure)

    newData <- new("metabData")

    newData = detectFields(Data = newData, table = table,
                           mz = mz, rt = rt, id, adduct = adduct,
                           samples = samples, extra = extra, Q = Q)

    newData = adjustData(Data = newData, misspc = misspc, measure = measure,
                         rtmin = rtmin, rtmax = rtmax, zero = zero,
                         impute = impute, imputeVal = imputeVal,
                         duplicate = duplicate)

    return(newData)
}






