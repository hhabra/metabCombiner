##
#' @title Filter Features by Retention Time Range
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


#' Duplicate Feature Detection Parameters
#'
#' @description Lists the parameters for detection of two or more rows that
#' represent the same entity, based on similar m/z and retention time values.
#'
#' @param mz m/z tolerance for duplicate feature detection
#'
#' @param rt RT tolerance for duplicate feature
#'
#' @param resolve character. Either "single" (default) or "merge".
#'
#' @param weighted logical. Option to weight m/z, RT, Q by mean abundance of
#' each row (TRUE) or take single representative values (FALSE).
#'
#' @details
#' The presence of duplicate features has negative consequences for the LC-MS
#' alignment task. The package offers several options for resolving the issue
#' of feature duplication. Pairwise m/z and RT tolerances define which features
#' are to be considered as duplicates within a single data set. Setting
#' \code{mz} or \code{rt} to 0 skips duplicate feature filtering altogether.
#'
#' When duplicates are detected, either a single master copy is retained
#' (\code{resolve} = "single") or merged into a single row
#' (\code{resolve} = "merge").The master copy is the copy with lower proportion
#' of missingness, followed by the most abundant (by median or mean). If %
#' missingness and abundance is equivalent for duplicates, the first copy that
#' appears is retained. The "merge" option fuses duplicate feature rows, with
#' quantitative descriptors (m/z, RT) either calculated as a weighted average
#' (\code{weighted} = TRUE) or otherwise taken from the top representative row;
#' id and adduct values are concatenated; the maximum feature value is used for
#' each sample; and all 'extra' values are taken from the 'master copy' row,
#' similar to the "single" option.
#'
#' @export
#'
#' @examples
#' data(plasma20)
#' pars.duplicate <- opts.duplicate(mz = 0.01, rt = 0.05, resolve = "single")
#' p20 <- metabData(plasma20, samples = "Red", duplicate = pars.duplicate)
#'
#' #to prevent removal of duplicate features
#' p20 <- metabData(plasma20, samples = "Red", duplicate = opts.duplicate(0))
#'
#' ##merge option
#' pars.duplicate <- opts.duplicate(mz = 0.01, rt = 0.05, resolve = "merge")
#' p20 <- metabData(plasma20, samples = "Red", duplicate = pars.duplicate)
#'
#'
opts.duplicate <- function(mz = 0.0025, rt = 0.05,
                           resolve = c("single", "merge"), weighted = FALSE)
{
    list(mz = mz, rt = rt, resolve = match.arg(resolve), weighted = weighted)
}


duplicate.pars <- function(duplicate)
{
    if(is.vector(duplicate) & length(duplicate) == 2){
        tolMZ <- duplicate[[1]]; tolRT <- duplicate[[2]]
        resolve <- "single"; weighted = FALSE
    }
    else{
        tolMZ <- duplicate$mz; tolRT <- duplicate$rt; resolve <- duplicate$resolve
        weighted <- duplicate$weighted
    }
    if(!is.numeric(tolMZ) | !is.numeric(tolRT) | tolMZ < 0 | tolRT < 0)
        stop("duplicate tolerances must be numeric and positive")

    return(list(tolMZ = tolMZ, tolRT = tolRT, resolve = resolve,
                weighted = weighted))
}

cat.strings <- function(x){
    x <- unique(x[!x %in% c(NA, "")])
    if(length(x) == 0) return("")
    else if(length(x) == 1) return(x)
    else
        paste0("{", paste(x, collapse = "; "), "}")
}

summaryCounts <- function(sampleData, measure)
{
    if(measure == "median")
        counts <- matrixStats::rowMedians(as.matrix(sampleData), na.rm = TRUE)
    else if(measure == "mean")
        counts <- matrixStats::rowMeans2(as.matrix(sampleData), na.rm = TRUE)
    counts <- ifelse(is.na(counts), 0, counts)
    return(counts)
}

sampleData <- function(Data, data, zero)
{
    samples <- data[getSamples(Data)]
    if(isTRUE(zero))  samples[samples == 0] <- NA
    return(samples)
}


groupDuplicates <- function(mzrt, tolMZ, tolRT)
{
    mzrt[["dupbin"]] <- .Call("binDuplicates",
                                mz = mzrt[["mz"]],
                                tolMZ = as.numeric(tolMZ),
                                PACKAGE = "metabCombiner")
    mzrt <- dplyr::arrange(mzrt, .data$dupbin, .data$missing,
                           desc(.data$counts))
    dupdata <- .Call("groupDuplicates",
                     mz = mzrt[["mz"]],
                     rt = mzrt[["rt"]],
                     tolMZ = as.numeric(tolMZ),
                     tolRT = as.numeric(tolRT),
                     dupbin = as.integer(mzrt[["dupbin"]]),
                     PACKAGE = "metabCombiner")
    mzrt[["dupgroup"]] = dupdata[[1]];
    mzrt[["dupremove"]] = dupdata[[2]];
    return(mzrt);
}

drop.dup.columns <- function(data)
{
    data[!names(data) %in% c("dupgroup", "dupremove")]
}

merge_duplicate_mzrtQ <- function(dup, nondup, weighted)
{
    if(weighted)
        mzrtQ_dup <- dup[c("dupgroup", "mz", "rt", "Q", "wt")] %>%
            dplyr::group_by(.data$dupgroup) %>%
            dplyr::summarize(mz = stats::weighted.mean(.data$mz, w = .data$wt),
                             rt = stats::weighted.mean(.data$rt, w = .data$wt),
                             Q = stats::weighted.mean(.data$Q, w = .data$wt))
    else
        mzrtQ_dup <- dup[c("dupgroup", "mz", "rt", "Q")][dup$dupremove == 0,]
    mzrtQ_nondup <- nondup[c("dupgroup","mz","rt","Q")]
    mzrtQ <- rbind.data.frame(mzrtQ_dup, mzrtQ_nondup) %>%
            dplyr::arrange(.data$dupgroup)
    return(mzrtQ)
}

merge_duplicate_idadd <-function(dup, nondup)
{
    idadd_dup <- stats::aggregate(x = dup[c("dupgroup", "id", "adduct")],
                            by = list(dup$dupgroup), FUN = cat.strings)
    idadd_nondup <- nondup[c("dupgroup","id","adduct")]
    idadd <- rbind.data.frame(idadd_dup[,-1], idadd_nondup) %>%
            dplyr::arrange(.data$dupgroup)
    return(idadd)
}

merge_duplicate_values <- function(dup, nondup, sampnames)
{
    n <- length(sampnames)
    ngroup <- length(unique(dup[["dupgroup"]]))
    values_dup <- dup[c("dupgroup", sampnames)]
    values_nondup <- nondup[c("dupgroup", sampnames)]
    values_vector <- as.vector(as.matrix(dup[sampnames]))
    dupgroup_vector <- rep(values_dup[["dupgroup"]], length(sampnames))
    merged_vector <- .Call("merge_duplicate_values",
                           values = as.numeric(values_vector),
                           dupgroup = dupgroup_vector,
                           ngroup = as.integer(ngroup),
                           n = n, PACKAGE = "metabCombiner")
    merged_values <- cbind.data.frame(unique(values_dup[["dupgroup"]]),
                        matrix(merged_vector, nrow = ngroup, ncol = n))
    colnames(merged_values) <- names(values_dup)
    merged_values <- rbind.data.frame(merged_values, values_nondup) %>%
                    dplyr::arrange(.data$dupgroup)
    return(merged_values[-1])
}

mergeDuplicates <- function(Data, data, samples, measure, zero, weighted)
{
    n <- ncol(samples)
    data <- data[with(data, order(`dupgroup`)),]
    data[["wt"]] <- matrixStats::rowSums2(as.matrix(samples), na.rm = TRUE) / n
    data[["wt"]] <- ifelse(data[["wt"]] > 0, data[["wt"]], 1e-5)
    rowID <- data[["rowID"]][!data$dupremove]
    groupsWithDups <- unique(data[["dupgroup"]][data[["dupremove"]] != 0])
    nondup <- data[!data[["dupgroup"]] %in% groupsWithDups,]
    dup <- data[data[["dupgroup"]] %in% groupsWithDups,]
    mzrtQ <- merge_duplicate_mzrtQ(dup, nondup, weighted)
    idadd <- merge_duplicate_idadd(dup, nondup)
    merged_values <- merge_duplicate_values(dup, nondup, names(samples))
    extra <- data[getExtra(Data)][!data$dupremove,]
    data <- data.frame(rowID = rowID,id = idadd$id,mz = mzrtQ$mz,rt = mzrtQ$rt,
                       adduct = idadd$adduct, Q = mzrtQ$Q, merged_values,
                       drop.dup.columns(extra), check.names = FALSE)
    samples <- sampleData(Data, data, zero)
    missing <- matrixStats::rowCounts(as.matrix(samples), value = NA)
    counts <- summaryCounts(samples, measure)
    return(list(data = data, counts = counts, missing = missing))
}

filterDuplicates <- function(Data, data, missing, counts)
{
    counts <- counts[data$dupremove == 0]
    missing <- missing[data$dupremove == 0]
    data <- drop.dup.columns(data[data$dupremove == 0,])
    return(list(data = data, counts = counts, missing = missing))
}


duplicateHandler <- function(Data, data, samples, filtered, duplicate,
                             measure, missing, counts, zero)
{
    pars <- duplicate.pars(duplicate)
    if(pars[["tolMZ"]] == 0 | pars[["tolRT"]] == 0 | nrow(data) == 0)
        return(list(data = data, missing = missing, counts = counts,
                    filtered = filtered))
    mzrt <- dplyr::select(data, .data$rowID,.data$mz,.data$rt) %>%
            dplyr::mutate(missing = as.numeric(missing),
                        counts = as.numeric(counts)) %>%
            dplyr::arrange(.data$mz)
    mzrt <- groupDuplicates(mzrt, pars[["tolMZ"]], pars[["tolRT"]])
    mzrt <- dplyr::arrange(mzrt, .data$rowID)
    data[c("dupgroup", "dupremove")] <- mzrt[c("dupgroup", "dupremove")]
    filtered[["duplicates"]] <- drop.dup.columns(data[data$dupremove == 1,])
    weighted <- pars[["weighted"]]
    if(pars[["resolve"]] == "merge" & sum(data[["dupremove"]]) > 0)
        dupData <- mergeDuplicates(Data, data, samples, measure, zero, weighted)
    else
        dupData <- filterDuplicates(Data, data, missing, counts)
    dupData[["filtered"]] <- filtered
    return(dupData)
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
#' 2) Removal of features with missingness percentage exceeding \code{misspc}
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
#' \code{\link{metabData}}, \code{\link{filterRT}}
#'
##
adjustData <- function(Data, misspc, measure, rtmin, rtmax, zero, duplicate)
{
    data <- getData(Data)
    stats <- list()
    filtered <- list()
    stats[["input_size"]] <- nrow(data)
    data <- filterRT(data, rtmin = rtmin, rtmax = rtmax)
    filtered[["rt"]] <- setdiff(getData(Data), data)
    stats[["filtered_by_rt"]] <- stats[["input_size"]] - nrow(data)
    samples <- sampleData(Data, data, zero)
    missing <- matrixStats::rowCounts(as.matrix(samples), value = NA)
    counts <- summaryCounts(samples, measure)
    dd <- duplicateHandler(Data, data, samples, filtered, duplicate, measure,
                           missing, counts, zero)
    filtered <- dd[["filtered"]]
    stats[["filtered_as_duplicates"]] <- nrow(filtered[["duplicates"]])
    keepIndices <- which((dd[["missing"]] / ncol(samples) * 100) <= misspc)
    stats[["filtered_by_missingness"]] <- nrow(dd$data) - length(keepIndices)
    filtered[["missing"]] <- dd$data[-keepIndices,]
    data <- dd$data[keepIndices,]
    counts <- dd$counts[keepIndices]
    if(nrow(data) == 0) stop("empty dataset following missingness filter")
    stats[["final_count"]] <- nrow(data)
    if(identical(data[["Q"]], rep(0, nrow(data))))
        data["Q"] <- round((rank(counts) - 0.5) / length(counts), 4)
    Data <- update_md(Data, data = data, stats = stats, filtered = filtered)
    return(Data)
}




