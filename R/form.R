findDuplicates <- function(data, missing, counts, duplicate)
{
    tolMZ <- duplicate[1];  tolRT <- duplicate[2]
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

#' @title Determine Averaging Options
#'
#' @param means logical vector or list of logicals
#'
#' @noRd
form_means <- function(means){
    if(is.list(means)){
        mean_mz <- means[["mz"]]; if(is.null(mean_mz)) mean_mz = FALSE
        mean_rt <- means[["rt"]]; if(is.null(mean_rt)) mean_rt = FALSE
        mean_Q <- means[["Q"]]; if(is.null(mean_Q)) mean_Q = FALSE
        means <- c(mean_mz, mean_rt, mean_Q)
    }
    else if(length(means) == 1)    means <- rep(means, 3)
    else if(length(means) == 3)    means <- means
    else
        stop("invalid value supplied to 'means' argument")

    return(means)
}

#' Reform Report with New Labels
#'
#' @description Helper function for labelRows(). After determining row
#' annotations, stitches together the metadata, the row annotations, and
#' the sample + extra value data.
#'
#' @param fields data.frame abridged combinedTable with metadata fields.
#'
#' @param values data.frame combinedTable with sample (+ extra) columns
#'
#' @param group0 data.frame consisting of "group 0" features (obtained using
#'  update_MC method) to be appended
#'
#' @noRd
formLabeledTable <- function(fields, values, remove, group0 = NULL)
{
    fields$labels <- ifelse(is.na(fields$labels), "", fields$labels)

    #option to eliminate rows labeled as removables
    if(remove == TRUE){
        keepRows <- which(fields[["labels"]] != "REMOVE")
        fields <- fields[keepRows,]
        values <- values[keepRows,]
    }

    #eliminating potential duplicate columns
    labnames <- c("labels", "subgroup", "alt", "resolveScore")
    if(any(labnames %in% utils::head(names(values)))){
        duplicates <- which(head(names(values)) %in% labnames)
        values <- values[,-duplicates]
    }
    cTable <- data.frame(fields, values, stringsAsFactors = FALSE,
                         check.names = FALSE)

    if(!is.null(group0)){
        if("resolveScore" %in% names(fields))
            group0$fields[["resolveScore"]] <- 0
        group0 <- cbind.data.frame(group0[["fields"]], group0[["values"]])
        cTable <- rbind.data.frame(cTable, group0)
    }

    return(cTable)
}


#' @title Form Single Dataset from Object
#'
#' @description When a metabCombiner object is supplied as input for a new
#' \code{metabCombiner} construction, information from one dataset or the mean
#' of all constituent datasets is chosen to represent that object. Only one row
#' is allowed for each feature, with duplicates discarded.
#'
#' @param object metabCombiner object
#'
#' @param data character dataset identifer
#'
#' @param means logical length 3 vector or list indicating whether to take means
#'              of m/z, RT, or Q, respectively;
#'
#' @param impute logical. Should features with missing values
#' be imputed (currently only mean-value imputation available)
#'
#' @details These are the steps executed in this function:
#'
#' 1) resolving conflicting feature pair alignment (FPA) subgroups
#'
#' 2) removal of all FPAs annotated as "REMOVE" in the labels column
#'
#' 3) Reforming the combinedTable and matching the corresponding featData
#'
#' 4) selection of meta-data from a single dataset
#'
#' @noRd
form_dataset <- function(object, data, means, rtOrder, impute)
{
    cTable <- combinedTable(object)
    fnames <- c(combinerNames(), "labels", "subgroup", "alt")
    if(length(setdiff(fnames, names(cTable))) > 0)
        stop("labelRows() has not been called on object")
    fields <- resolveRows(cTable[fnames], rtOrder = rtOrder)
    values <- cTable[-seq(1,length(fnames))]
    cTable <- formLabeledTable(fields, values, remove = TRUE, group0 = NULL)
    fulldata <- featData(object)
    fulldata <- fulldata[match(cTable[["rowID"]], fulldata[["rowID"]]),]
    fdata <- fulldata[grep(paste("rowID|_", data, "$", sep = ""),
                           names(fulldata), value = TRUE)]
    names(fdata) <- c("rowID","id","mz","rt","Q","adduct")
    samps_extras <- rmbrackets(paste(unlist(lapply(datasets(object), function(d)
        c(getSamples(object, d), getExtra(object,d)))), collapse = "|"))
    vnames <- names(values)[grep(samps_extras, rmbrackets(names(values)))]
    dataset <- cbind.data.frame(fdata, cTable[,names(cTable) %in% vnames])
    ave_mz <- matrixStats::rowMeans2(as.matrix(
        fulldata[grep("^mz_", names(fulldata), value = TRUE)]), na.rm = TRUE)
    ave_rt <- matrixStats::rowMeans2(as.matrix(
        fulldata[grep("^rt_", names(fulldata), value = TRUE)]), na.rm = TRUE)
    ave_Q <- matrixStats::rowMeans2(as.matrix(
        fulldata[grep("^Q_", names(fulldata), value = TRUE)]), na.rm = TRUE)
    naInds <- which(is.na(dataset$mz) | is.na(dataset$rt) | is.na(dataset$Q))
    if(length(naInds) > 0){
        if(isTRUE(impute)){
            dataset$mz <- ifelse(is.na(dataset$mz), ave_mz, dataset$mz)
            dataset$rt <- ifelse(is.na(dataset$rt), ave_rt, dataset$rt)
            dataset$Q <- ifelse(is.na(dataset$Q), ave_Q, dataset$Q)
        }
        else{
            dataset <- dataset[-naInds,]
            ave_mz <- ave_mz[-naInds]
            ave_rt <- ave_rt[-naInds]
            ave_Q <- ave_Q[-naInds]
        }
    }
    means <- form_means(means)
    if(isTRUE(means[1]))  fdata[["mz"]] <- ave_mz
    if(isTRUE(means[2]))  fdata[["rt"]] <- ave_rt
    if(isTRUE(means[3]))  fdata[["Q"]] <- ave_Q
    duplicates <- findDuplicates(dataset, missing = rep(0, nrow(dataset)),
                                counts = dataset$Q, duplicate = c(1E-4, 1E-4))
    if(length(duplicates) > 0) dataset <- dataset[-duplicates,]
    return(dataset)
}

#' @title Form metabCombiner Report Table
#'
#' @description  Takes previously computed m/z groups using \code{mzGroup()}
#' and creates merged \code{combinedTable} consisting of all possible feature
#' alignments, all initialized with equivalent score and ranking.
#'
#' @param object metabCombiner object
#'
#' @param xset data frame. A processed metabolomics feature table.
#'
#' @param yset data frame. A processed metabolomics feature table.
#'
#' @param xreps number of repeats for xset
#'
#' @param yreps number of repeats for yset
#'
#' @return metabCombiner object with initialized combinedTable data frame.
#'
#' @noRd
##
formCombinedTable <- function(object, xset, yset, xreps, yreps){
    xComb <- dplyr::slice(xset, rep(seq(1,n()), times = xreps))
    yComb <- yset[yreps,]
    samps_extras <- unlist(lapply(datasets(object), function(d)
        c(getSamples(object, d), getExtra(object,d))))
    x_samps_extras <- unique(stats::na.omit(match(samps_extras, names(xComb))))
    x_samps_extras <- x_samps_extras[x_samps_extras %in% seq(7, ncol(xComb)-1)]
    y_samps_extras <- unique(stats::na.omit(match(samps_extras, names(yComb))))
    y_samps_extras <- y_samps_extras[y_samps_extras %in% seq(7, ncol(yComb)-1)]
    cTable <- data.frame(idx = xComb[["id"]], idy = yComb[["id"]],
                    mzx = round(xComb[["mz"]],5), mzy = round(yComb[["mz"]],5),
                    rtx = round(xComb[["rt"]],4), rty = round(yComb[["rt"]],4),
                    rtProj = numeric(nrow(xComb)), Qx = round(xComb[["Q"]],4),
                    Qy = round(yComb[["Q"]],4), group = xComb[["group"]],
                    score = rep(1, nrow(xComb)), rankX = as.integer(1),
                    rankY = as.integer(1), adductx = xComb[["adduct"]],
                    adducty = yComb[["adduct"]], xComb[,sort(x_samps_extras)],
                    yComb[,sort(y_samps_extras)], stringsAsFactors = FALSE,
                    check.names = FALSE)
    consec <- lapply(rle(cTable[["group"]])[["lengths"]], function(l) seq(1,l))
    consec <- unlist(consec)
    rowID <- sprintf("%s.%s", cTable[["group"]], consec)
    cTable <- cbind.data.frame(rowID = rowID, cTable)

    return(cTable)
}

#' @title Form Feature Meta-data
#'
#' @description This forms the featData slot of a metabCombiner object,
#' consisting of feature descriptors from all constituent datasets.
#'
#' @param object metabCombiner object
#'
#' @param xfeat data frame of X dataset feature descriptors
#'
#' @param yfeat data frame of Y dataset feature descriptors
#'
#' @param xreps Number of repeats for X dataset features
#'
#' @param yreps Number of repeat for Y dataset features
#'
#' @noRd
form_featData <- function(object, cTable, xfeat, yfeat, xreps, yreps)
{
    if(is.null(xfeat) & is.null(yfeat))
        featData <- cbind.data.frame(cTable[["idx"]], cTable[["idy"]],
                                     cTable[["mzx"]], cTable[["mzy"]],
                                     cTable[["rtx"]], cTable[["rty"]],
                                     cTable[["Qx"]], cTable[["Qy"]],
                                     cTable[["adductx"]], cTable[["adducty"]])
    else if(!is.null(xfeat) & is.null(yfeat)){
        xfeat <- dplyr::slice(xfeat, rep(seq(1,n()), times = xreps))
        featData <- cbind.data.frame(
            dplyr::select(xfeat, starts_with("id_")), cTable[["idy"]],
            dplyr::select(xfeat, starts_with("mz_")), cTable[["mzy"]],
            dplyr::select(xfeat, starts_with("rt_")), cTable[["rty"]],
            dplyr::select(xfeat, starts_with("Q_")),  cTable[["Qy"]],
            dplyr::select(xfeat, starts_with("adduct_")),cTable[["adducty"]])
    }
    else if(is.null(xfeat) & !is.null(yfeat)){
        yfeat <- yfeat[yreps,]
        featData <- cbind.data.frame(
            cTable[["idx"]], dplyr::select(yfeat, starts_with("id_")),
            cTable[["mzx"]], dplyr::select(yfeat, starts_with("mz_")),
            cTable[["rtx"]], dplyr::select(yfeat, starts_with("rt_")),
            cTable[["Qx"]],  dplyr::select(yfeat, starts_with("Q_")),
            cTable[["adductx"]],dplyr::select(yfeat, starts_with("adduct_")))
    }
    else if(!is.null(xfeat) & !is.null(yfeat)){
        xfeat <- dplyr::slice(xfeat, rep(seq(1,n()), times = xreps))
        yfeat <- yfeat[yreps,]
        featData <- cbind.data.frame(
                                dplyr::select(xfeat, starts_with("id_")),
                                dplyr::select(yfeat, starts_with("id_")),
                                dplyr::select(xfeat, starts_with("mz_")),
                                dplyr::select(yfeat, starts_with("mz_")),
                                dplyr::select(xfeat, starts_with("rt_")),
                                dplyr::select(yfeat, starts_with("rt_")),
                                dplyr::select(xfeat, starts_with("Q_")),
                                dplyr::select(yfeat, starts_with("Q_")),
                                dplyr::select(xfeat, starts_with("adduct_")),
                                dplyr::select(yfeat, starts_with("adduct_")))
    }
    fields <- c("id","mz","rt","Q","adduct")
    names(featData) <- paste(rep(fields, each = length(datasets(object))),
                             rep(datasets(object),length(fields)), sep = "_")
    featData <- cbind.data.frame(rowID = cTable[["rowID"]], featData)
    return(featData)
}

#' @title Form combinedTable and featData
#'
#' @description This wraps together the formation of the combinedTable and
#' featData from the input datasets used to construct the metabCombiner object.
#' This completes the initial object formation from the features previously
#' grouped by m/z.
#'
#' @param object metabCombiner object
#'
#' @param xset data frame of grouped X dataset features used to form the
#'              combinedTable result
#'
#' @param yset data frame of grouped X dataset features used to form the
#'              combinedTable result
#'
#' @param xfeat data frame of feature metadata originating from xdata input
#'
#' @param yfeat data frame of feature metadata originating from ydata input
#'
#' @param nGroups total number of feature groups as determined by similar m/z
#'
#'
#' @noRd
form_tables <- function(object, xset, yset, xfeat = NULL, yfeat = NULL, nGroups)
{
    groupCountX <- as.integer(table(xset[["group"]]))
    groupCountY <- as.integer(table(yset[["group"]]))
    xreps <- rep(groupCountY, times = groupCountX)
    ends <- c(which(!duplicated(yset[["group"]])),nrow(yset)+1)
    yreps <- unlist(lapply(seq(1,nGroups), function(number){
        rep(seq(ends[number], ends[number+1]-1), times = groupCountX[number])
    }))
    cTable <- formCombinedTable(object, xset, yset, xreps, yreps)
    featData <- form_featData(object, cTable, xfeat, yfeat, xreps, yreps)
    object <- update_mc(object, combinedTable = cTable, featData = featData)
    return(object)
}
