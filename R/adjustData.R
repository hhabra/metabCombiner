##
#' @title Filter Features by Retention Time
#'
#' @description
#' Restricts input metabolomics feature table in \code{metabData} object to a
#' range of retention times defined by \code{rtmin} & \code{rtmax}.
#'
#' @param data   formatted metabolomics data frame.
#'
#' @param rtmin  lower range of retention times for analysis. If "min", defaults
#' to minimum observed retention time.
#' .
#' @param rtmax  upper range of retention times for analysis. If "max", defaults
#' to maximum observed retention time.
#'
#' @details
#' Retention time restriction is often a recommended step to aid the analysis of
#' comparable metabolomics datasets. For example, the beginning and end of a
#' chromatogram often contain features that do not correspond with a true
#' a biological compound derived from the sample. It is recommended to set rtmin
#' slightly before the first observed common metabolite, and rtmax slightly after
#' the last observed metabolite.
#'
#' @return A data frame of metabolomics features, limited to retention time window
#' \code{rtmin} \eqn{\le} rt \eqn{\le} \code{rtmax})
#'
##
filterRT <- function(data, rtmin, rtmax)
{
    rts = data[["rt"]]

    if(rtmin == "min")
        rtmin = min(rts)

    if(rtmax == "max")
        rtmax = max(rts)

    if(class(rtmin)!= "numeric" | rtmin > rtmax | rtmin < 0 |
        rtmin > max(rts) | rtmin < min(rts))
    {
        warning("The supplied 'rtmin value' is invalid. Setting 'rtmin' to
                minimum observed retention time.")
        rtmin = min(rts)
    }

    if(class(rtmax)!= "numeric" | rtmin > rtmax | rtmax < 0 | rtmax > max(rts)){
        warning("The supplied 'rtmax' value is invalid. Setting to maximum
                 observed retention time.")
        rtmax = max(rts)
    }

    data = dplyr::filter(data, rt >= rtmin & rt <= rtmax)

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
#' Pairs of features are deemed duplicates if their differences in m/z & rt fall
#' within tolerances defined by the \code{duplicate} parameter. If a pair of
#' duplicate features is found, one of the pair is removed. The determination of
#' which feature to remove is first by percent missingness, followed by which has
#' a lower central abundance measure. If the features have equal missingness and
#' abundance, then row priority arbitrarily determines the feature to be removed.
#'
#' @return  integer indices of removable duplicate features
##
findDuplicates <- function(data, missing, counts, duplicate)
{
    if(length(duplicate)!= 2)
        stop("Argument 'duplicate' must be a numeric, positive ordered pair")

    tolMZ = duplicate[1]
    tolRT = duplicate[2]

    if(!is.numeric(tolMZ) | !is.numeric(tolRT))
        stop("Parameter 'duplicate' must be a numeric, positive ordered pair")

    if(tolMZ < 0 | tolRT < 0)
        stop("Parameter 'duplicate' must be a numeric, positive ordered pair")

    if(tolMZ == 0 | tolRT == 0)
        return(numeric(0))

    dataMat = dplyr::select(data, mz,rt) %>%
              dplyr::mutate(counts = counts,
                            missing = missing,
                            index = 1:nrow(data)) %>%
              dplyr::arrange(mz)

    dataMat[["labels"]] = .Call("findDuplicates",
                                mz = dataMat[["mz"]],
                                rt = dataMat[["rt"]],
                                tolMZ = tolMZ,
                                tolRT = tolRT,
                                missing = dataMat[["missing"]],
                                counts = dataMat[["counts"]],
                                PACKAGE = "Combiner")

    duplicates = dataMat[["index"]][dataMat[["labels"]] == 1]

    return(duplicates)
}

##
#' Process and Filter Metabolomics Feature Lists
#'
#' @description
#' \code{adjustData()} contains a set of pre-analysis steps for processing LC-MS
#' metabolomics feature tables.
#'
#' @param Data   a metabData object.
#'
#' @param misspc  Numeric. Threshold missingness percentage for analysis.
#'
#' @param measure  Character. Central abundance measure, either "median" or "mean".
#'
#' @param rtmin   Numeric. Minimum retention time for analysis.
#'
#' @param rtmax   Numeric. Maximum retention time for analysis.
#'
#' @param zero  Logical. Whether to consider zero values as missing.
#'
#' @param duplicate Ordered numeric pair (m/z, rt) tolerance parameters for
#' duplicate feature search.
#'
#' @details
#' The pre-analysis adjustment steps include: 1) Restriction to a feature
#' retention time range \code{rtmin} \eqn{\le} rt \eqn{\le} \code{rtmax}
#' 2) Removal of features with a missing value percentage exceeding \code{misspc}
#' 3) Removal of duplicate metabolomics features.
#'
#' After reducing the table, abundance quantile values are calculated for the
#' remaining features (unless initially supplied by the user). The \code{measure}
#' argument determines which central abundance statistic ("mean" or "median") is
#' used to rank features.
#'
#' @return
#' An updated metabData object. The \code{data} field is processed by the listed
#' steps and
#'
#' @seealso
#' \code{\link{metabData}}: the constructor for metabData objects,
#' \code{\link{filterRT}}: function for filtering by retention times,
#' \code{\link{findDuplicates}}: function for finding duplicates
#'
##
adjustData <- function(Data, misspc, measure, rtmin, rtmax, zero,duplicate)
{
    data = getData(Data)
    samples = getSamples(Data)

    stats = list()
    stats[["input_size"]] = nrow(data)

    ##filtering by RT region
    data = filterRT(data, rtmin = rtmin, rtmax = rtmax)

    stats[["filtered_by_rt"]] = stats[["input_size"]] - nrow(data)

    ##filtering by % missingness
    missingpc <- apply(data[samples], 1, function(row){
        if(zero == TRUE)
            row[row <= 0] <- NA

        sum(is.na(row)) / length(row) * 100
    })

    keepIndices = which(missingpc <= misspc)

    stats[["filtered_by_missingness"]] = nrow(data) - length(keepIndices)

    data = data[keepIndices,]
    missingpc = missingpc[keepIndices]

    if(measure == "median")
        counts <- apply(data[,samples], 1, stats::median, na.rm = TRUE) %>%
                  as.numeric()
    else if(measure == "mean")
        counts <- apply(data[,samples], 1, mean, na.rm = TRUE) %>%
                  as.numeric()

    ##removing duplicate features
    duplicates = findDuplicates(data = data,
                                counts = counts,
                                missing = missingpc,
                                duplicate = duplicate)

    if(length(duplicates) > 0){
        data = data[-duplicates,]
        counts = counts[-duplicates]
    }

    ##optional imputation of missing values
    if(impute == TRUE & any(missingpc > 0)){
        data = imputeVals(data = data,
                          samples = samples,
                          imputeVal = imputeVal)
    }

    stats[["filtered_as_duplicates"]] = length(duplicates)
    stats[["final_count"]] = nrow(data)

    ##calculating abundance quantiles
    if(identical(data[["Q"]], rep(0, nrow(data))))
        data["Q"] <- round((rank(counts) - 0.5) / length(counts),4)

    Data@data = data
    Data@stats = stats

    return(Data)
}
