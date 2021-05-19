##
#' @title Filter Features by Retention Time
#'
#' @description
#' Restricts input metabolomics feature table in \code{metabData} object to a
#' range of retention times defined by \code{rtmin} & \code{rtmax}.
#'
#' @param data   formatted metabolomics data frame.
#'
#' @param rtmin  lower range of retention times for analysis. If "min",
#' defaults to minimum observed retention time.
#' .
#' @param rtmax upper range of retention times for analysis. If "max", defaults
#' to maximum observed retention time.
#'
#' @details
#' Retention time restriction is often recommended to aid the analysis of
#' comparable metabolomics datasets. The beginning and end of a chromatogram
#' typically contain features that do not correspond with true biological
#' compounds derived from the sample. \code{rtmin} and \code{rtmax} should be
#' set slightly before and slightly after the first and last commonly observed
#' metabolites, respectively.
#'
#'
#' @return A data frame of metabolomics features, limited to time window
#' \code{rtmin} \eqn{\le} rt \eqn{\le} \code{rtmax})
#'
##
filterRT <- function(data, rtmin, rtmax)
{
    rts <- data[["rt"]]

    if(rtmin == "min")
        rtmin <- min(rts)
    if(rtmax == "max")
        rtmax <- max(rts)
    if(!is.numeric(rtmin) | rtmin > rtmax | rtmin < 0 | rtmin > max(rts))
        stop("invalid value supplied for 'rtmin' argument")
    if(!is.numeric(rtmax) | rtmax < min(rts) | rtmax < 0)
        stop("invalid value supplied for 'rtmax' argument")
    data <- dplyr::filter(data, .data$rt >= rtmin & .data$rt <= rtmax)

    if(nrow(data) == 0)
        stop("empty dataset following retention time filter")

    return(data)
}

##
#' @title Find and Remove Duplicate Features
#'
#' @description
#' Pairs of features with nearly identical m/z and retention time values are
#' removed in this step.
#'
#' @param data     Constructed metabolomics data frame.
#'
#' @param duplicate   Ordered numeric pair (m/z, rt) tolerance parameters for
#' duplicate feature search.
#'
#' @param missing  Numeric vector. Percent missingness for each feature.
#'
#' @param counts   Numeric vector. Central measure for each feature.
#'
#' @details
#' Pairs of features are deemed duplicates if pairwise differences in m/z & rt
#' fall within tolerances defined by the \code{duplicate} argument. If a pair
#' of duplicate features is found, one member is removed. The determination of
#' which feature to remove is first by percent missingness, followed by central
#' abundance measure (median or mean). If the features have equal missingness
#' and abundance, then row order determines the feature to be removed.
#'
#' @return  integer indices of removable duplicate features
##
findDuplicates <- function(data, missing, counts, duplicate)
{
    if(length(duplicate)!= 2)
        stop("'duplicate' argument must be a numeric, positive ordered pair")

    tolMZ <- duplicate[1]
    tolRT <- duplicate[2]

    if(!is.numeric(tolMZ) | !is.numeric(tolRT) | tolMZ < 0 | tolRT < 0)
        stop("'duplicate' argument must be a numeric, positive ordered pair")

    if(tolMZ == 0 | tolRT == 0 | nrow(data) == 0)
        return(numeric(0))

    datMatrix <- dplyr::select(data, .data$mz,.data$rt) %>%
                dplyr::mutate(counts = counts,
                            missing = missing,
                            index = seq(1,nrow(data))) %>%
                dplyr::arrange(.data$mz)

    datMatrix[["labels"]] <- .Call("findDuplicates",
                                mz = datMatrix[["mz"]],
                                rt = datMatrix[["rt"]],
                                tolMZ = as.numeric(tolMZ),
                                tolRT = as.numeric(tolRT),
                                missing = as.numeric(datMatrix[["missing"]]),
                                counts = as.numeric(datMatrix[["counts"]]),
                                PACKAGE = "metabCombiner")

    duplicates <- datMatrix[["index"]][datMatrix[["labels"]] == 1]

    return(duplicates)
}

##
#' Process and Filter Metabolomics Feature Lists
#'
#' @description
#' \code{adjustData} contains a set of pre-analysis steps for processing LC-MS
#' metabolomics feature tables individually
#'
#' @param Data a metabData object.
#'
#' @param misspc Numeric. Threshold missingness percentage for analysis.
#'
#' @param measure Character. Choice of central abundance measure; either
#' "median" or "mean".
#'
#' @param rtmin Numeric. Minimum retention time for analysis.
#'
#' @param rtmax Numeric. Maximum retention time for analysis.
#'
#' @param zero Logical. Whether to consider zero values as missing.
#'
#' @param duplicate Ordered numeric pair (m/z, rt) tolerance parameters for
#' duplicate feature search.
#'
#' @details
#' The pre-analysis adjustment steps include:
#' 1) Restriction to a feature retention time range \code{rtmin} \eqn{\le}
#'    rt \eqn{\le} \code{rtmax}
#' 2) Removal of features with a percent missingness exceeding \code{misspc}
#' 3) Removal of duplicate metabolomics features.
#'
#' After processing, abundance quantile (Q) values are calculated between 0 & 1
#' for the remaining features, as ranked by the \code{measure} argument, unless
#' provided by the user.
#'
#' @return
#' Updated \code{metabData} object. The \code{data} field is processed by the
#' listed steps and \code{stats} list updated to contain feature statistics.
#'
#' @seealso
#' \code{\link{metabData}}: the constructor for metabData objects,
#' \code{\link{filterRT}}: function for filtering by retention times,
#' \code{\link{findDuplicates}}: function for finding duplicates
#'
##
adjustData <- function(Data, misspc, measure, rtmin, rtmax, zero,duplicate)
{
    data <- getData(Data)
    samples <- getSamples(Data)
    stats <- list()
    stats[["input_size"]] <- nrow(data)

    data <- filterRT(data, rtmin = rtmin, rtmax = rtmax)
    stats[["filtered_by_rt"]] = stats[["input_size"]] - nrow(data)

    mpc <- matrixStats::rowCounts(as.matrix(data[samples]), value = NA)
    if(zero == TRUE)
        mpc <- mpc + matrixStats::rowCounts(as.matrix(data[samples]),
                                            value = 0, na.rm = TRUE)
    keepIndices <- which((mpc / length(samples) * 100) <= misspc)
    stats[["filtered_by_missingness"]] <- nrow(data) - length(keepIndices)
    data <- data[keepIndices,]
    if(nrow(data) == 0)
        stop("empty dataset following missingness filter")
    mpc <- mpc[keepIndices] / length(samples) * 100

    if(measure == "median")
        counts <- matrixStats::rowMedians(as.matrix(data[samples]),na.rm = TRUE)
    else if(measure == "mean")
        counts <- matrixStats::rowMeans2(as.matrix(data[samples]),na.rm = TRUE)

    duplicates <- findDuplicates(data = data, counts = counts, missing = mpc,
                                duplicate = duplicate)
    if(length(duplicates) > 0){
        data <- data[-duplicates,]
        counts <- counts[-duplicates]
    }
    stats[["filtered_as_duplicates"]] <- length(duplicates)
    stats[["final_count"]] <- nrow(data)

    if(identical(data[["Q"]], rep(0, nrow(data))))
        data["Q"] <- round((rank(counts) - 0.5) / length(counts),4)

    Data <- update_md(Data, data = data, stats = stats)
    return(Data)
}
