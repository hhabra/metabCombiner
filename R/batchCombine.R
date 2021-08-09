#' @title Stepwise Multi-batch LC-MS Alignment
#'
#' @description This is a method for aligning multiple batches of a single
#' metabolomics experiment in a stepwise manner using the \code{metabCombiner}
#' workflow. The input is a list of \code{metabData} objects corresponding to
#' the batch data frames arranged in sequential order (i.e. batch 1,2,...,N),
#' and parameter lists for each step; the output is an aligned feature table and
#' a \code{metabCombiner} object composed from the input batches.
#'
#' @param batches list of metabData objects corresponding to each LC-MS batch
#'
#' @param binGap numeric parameter used for grouping features by m/z.
#' See ?mzGroup for more details.
#'
#' @param fitMethod RT spline-fitting method, either "gam" or "loess"
#'
#' @param means logical. Option to take average m/z, rt, and/or Q from
#' \code{metabComber}. May be a 3-length vector, single value (TRUE/FALSE),
#' or a list with names "mz", "rt", "Q" as names.
#'
#' @param anchorParam list of parameter values for selectAnchors() function
#'
#' @param fitParam list of parameter values for fit_gam() or fit_loess()
#'
#' @param scoreParam list of parameter values for calcScores()
#'
#' @param reduceParam list of parameter values for reduceTable()
#'
#' @details
#' Retention time drifting is commonly observed in large-scale LC-MS experiments
#' in which samples are analyzed in multiple batches. Conventional LC-MS
#' pre-processing approaches may effectively align features detected in samples
#' from within a single batch, but fail in many cases to account for inter-batch
#' drifting, leading to misaligned features.
#'
#' \code{batchCombine} assumes that each batch has been previously processed
#' separately using conventional LC-MS preprocessing approaches (e.g. XCMS),
#' and can be represented as a data frame. Each batch data feature table must be
#' filtered and formatted as a \code{metabData} object and the batches must be
#' arranged as a list in sequential order of acquisition.
#'
#' \code{batchCombine} applies the \code{metabCombine} wrapper function to
#' successive pairs of metabolomics batches in a stepwise manner. Each iteration
#' consists of the key steps in the package workflow (feature m/z grouping,
#' anchor selection, retention time spline fitting, pairwise scoring, & table
#' reduction). The first two batches are aligned together, then the combined
#' results are aligned with the third batch, and so forth. Parameters for each
#' sub-method are arranged in list format, with their respective defaults (e.g.
#' fitgamParam() lists the default values for the fit_gam function).
#'
#' Following each iteration, m/z, rt, and Q values from the combined dataset may
#' be averaged to use for comparison with the next batch's feature quantitative
#' descriptors, if the \code{means} argument is set to TRUE; if set to FALSE,
#' feature information is drawn from the latter of the  previously combined
#' batches, identical to the manner in which id & adduct descriptors are drawn.
#'
#' @return
#' \item{object}{metabCombiner object of the final alignment; x is set to the
#'               penultimate batch and y is set to the final batch}
#' \item{table}{combined feature table consisting of feature descriptor values
#'              followed by per-sample abundances and extra columns}
#'
#' @note \code{batchCombine} is designed for aligning multi-batch datasets, i.e.
#' where each batch is acquired in a roughly identical manner. It is not
#' for disparately acquired LC-MS datasets (e.g. from different instruments,
#' chromatographic systems, laboratories, etc.).
#'
#' @examples
#' \donttest{
#' #identically formatted batches in list form
#' data(metabBatches)
#'
#' #obtain list of metabData objects
#' batchdata <- lapply(metabBatches, metabData, samples = "POOL",
#'                     extra = "SAMP", zero = TRUE)
#'
#' #optional: give each batch dataset a name
#' names(batchdata) <- paste("B", seq_along(batchdata), sep = "")
#'
#' #customize main workflow parameter lists
#' saparam <- selectAnchorsParam(tolmz = 0.002, tolQ = 0.2, tolrtq = 0.1)
#' fgparam <- fitgamParam(k = 20, iterFilter = 1)
#' csparam <- calcScoresParam(A = 70, B = 35, C = 0.3)
#' rtparam <- reduceTableParam(minScore = 0.5, maxRTerr = 0.33)
#'
#' #run batchCombine program
#' combinedRes <- batchCombine(batches = batchdata, binGap = 0.0075,
#'                means = list('mz' = TRUE, 'rt' = FALSE, 'Q' = FALSE),
#'                anchorParam = saparam, fitParam = fgparam,
#'                scoreParam = csparam, reduceParam = rtparam)
#'
#' #aligned table results & metabCombiner object results
#' cTable <- combinedRes$table
#' object <- combinedRes$object
#'
#' #if names were set earlier, the names should be returned by this
#' datasets(object)
#'}
#'
#' @seealso \code{\link{metabCombine}}
#'
#' @export
batchCombine <- function(batches, binGap = 0.005, fitMethod = "gam",
                        means = list('mz' = TRUE, 'rt' = TRUE, 'Q' = TRUE),
                        anchorParam = selectAnchorsParam(),
                        fitParam = fitgamParam(),
                        scoreParam = calcScoresParam(B = 30),
                        reduceParam = reduceTableParam())
{
    ts <- vapply(batches, function(d) methods::is(d, "metabData"), logical(1))
    if(any(ts == FALSE))
        stop("'batches' must be a list consisting of metabData objects")
    if(is.null(names(batches)))
        names(batches) <- as.character(seq(1, length(batches)))
    if(any(duplicated(names(batches))))
        stop("duplicate names detected in batches argument")
    reduceParam[["remove"]] <- TRUE;   reduceParam[["method"]] <- "score"
    reduceParam[["resolveConflicts"]] <- TRUE
    object <- batches[[1]]
    for(b in seq(2,length(batches))){
        xid <- names(batches)[b-1];  yid <- names(batches)[b]
        message("Aligning ", xid, " with ", yid)
        object <- metabCombine(xdata = object, ydata = batches[[b]],
                            xid = xid, yid = yid, means = means,
                            fitMethod = fitMethod, fitParam = fitParam,
                            anchorParam =  anchorParam,
                            scoreParam = scoreParam, labelParam = reduceParam)
    }
    samples_extras <- unlist(lapply(datasets(object), function(d)
        c(getSamples(object, d), getExtra(object, d))))
    cTable <- combinedTable(object)[samples_extras]
    table <- cbind.data.frame(featdata(object), cTable)
    output <- list(object = object, table = table)

    return(output)
}




