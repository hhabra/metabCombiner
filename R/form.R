#'
#'
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
#' @noRd
formLabeledTable <- function(fields, values, remove)
{
    #option to eliminate rows labeled as removables
    if(remove == TRUE){
        keepRows <- which(fields[["labels"]] != "REMOVE")
        fields <- fields[keepRows,]
        values <- values[keepRows,]
    }

    #eliminating potential duplicate columns
    labnames <- c("labels", "subgroup", "alt", "resolveScore")
    if(any(labnames %in% utils::head(names(values)))){
        duplicates <- grep(paste(labnames, collapse = "|"), head(names(values)))
        values <- values[,-duplicates]
    }

    cTable <- data.frame(fields, values, stringsAsFactors = FALSE,
                         check.names = FALSE)

    return(cTable)
}


#' @title Form Single Dataset from Object
#'
#' @description When a metabCombiner object is supplied as input for a new
#' metabCombiner construction, information from one dataset or the mean of all
#' constituent datasets is chosen to represent that object. Only one row is
#' allowed for each feature, with duplicates discarded.
#'
#' @param object metabCombiner object
#'
#' @param data character dataset identifer
#'
#' @param means logical length 3 vector or list indicating whether to take means
#'              of m/z, RT, or Q, respectively;
#'
#' @noRd
form_dataset <- function(object, data, means, rtOrder)
{
    cTable <- combinedTable(object)
    fields <- resolveRows(cTable[,seq(1,19)], rtOrder)
    values <- cTable[,seq(20,ncol(cTable))]
    cTable <- formLabeledTable(fields, values, remove = TRUE)

    fdata <- featdata(object, data = data)
    fdata <- fdata[match(cTable[["rowID"]], fdata[["rowID"]]),]

    samples_extras <- unlist(lapply(datasets(object), function(d)
        c(getSamples(object, d), getExtra(object,d))))
    names(fdata) <- c("rowID", "id","mz","rt","Q","adduct")
    means <- form_means(means)
    if(means[1] == TRUE)  fdata[["mz"]] <- mzdata(object, value = "mean")
    if(means[2] == TRUE)  fdata[["rt"]] <- rtdata(object, value = "mean")
    if(means[3] == TRUE)  fdata[["Q"]] <- Qdata(object, value = "mean")
    dataset <- cbind.data.frame(fdata, cTable[samples_extras])
    duplicates <- findDuplicates(dataset, missing = rep(0, nrow(dataset)),
                                 counts = dataset$Q, duplicate = c(1E-4, 1E-4))
    if(length(duplicates) > 0)
        dataset <- dataset[-duplicates,]
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
    yComb = yset[yreps,]
    samples_extras <- unlist(lapply(datasets(object), function(d)
        c(getSamples(object, d), getExtra(object,d))))
    x_samples_extras <- intersect(names(xComb), samples_extras)
    y_samples_extras <- intersect(names(yComb), samples_extras)
    cTable <- data.frame(idx = xComb[["id"]], idy = yComb[["id"]],
                    mzx = round(xComb[["mz"]],5), mzy = round(yComb[["mz"]],5),
                    rtx = round(xComb[["rt"]],4), rty = round(yComb[["rt"]],4),
                    rtProj = numeric(nrow(xComb)),
                    Qx = round(xComb[["Q"]],4), Qy = round(yComb[["Q"]],4),
                    group = xComb[["group"]],
                    score = rep(1, nrow(xComb)),
                    rankX = as.integer(1), rankY = as.integer(1),
                    adductx = xComb[["adduct"]], adducty = yComb[["adduct"]],
                    xComb[x_samples_extras], yComb[y_samples_extras],
                    stringsAsFactors = FALSE, check.names = FALSE)
    consec <- lapply(rle(cTable[["group"]])[["lengths"]], function(l) seq(1,l))
    consec <- unlist(consec)
    rowID <- sprintf("%s.%s", cTable[["group"]], consec)
    cTable <- cbind.data.frame(rowID = rowID, cTable)

    return(cTable)
}

#' @title Form Feature Metadata
#'
#' @description This forms the featdata slot of a metabCombiner object,
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
form_featdata <- function(object, cTable, xfeat, yfeat, xreps, yreps)
{
    if(is.null(xfeat) & is.null(yfeat))
        featdata <- cbind.data.frame(cTable[["idx"]], cTable[["idy"]],
                                     cTable[["mzx"]], cTable[["mzy"]],
                                     cTable[["rtx"]], cTable[["rty"]],
                                     cTable[["Qx"]], cTable[["Qy"]],
                                     cTable[["adductx"]], cTable[["adducty"]])
    else if(!is.null(xfeat) & is.null(yfeat)){
        xfeat <- dplyr::slice(xfeat, rep(seq(1,n()), times = xreps))
        featdata <- cbind.data.frame(
            dplyr::select(xfeat, starts_with("id_")), cTable[["idy"]],
            dplyr::select(xfeat, starts_with("mz_")), cTable[["mzy"]],
            dplyr::select(xfeat, starts_with("rt_")), cTable[["rty"]],
            dplyr::select(xfeat, starts_with("Q_")),  cTable[["Qy"]],
            dplyr::select(xfeat, starts_with("adduct_")),cTable[["adducty"]])
    }
    else if(is.null(xfeat) & !is.null(yfeat)){
        yfeat <- yfeat[yreps,]
        featdata <- cbind.data.frame(
            cTable[["idx"]], dplyr::select(yfeat, starts_with("id_")),
            cTable[["mzx"]], dplyr::select(yfeat, starts_with("mz_")),
            cTable[["rtx"]], dplyr::select(yfeat, starts_with("rt_")),
            cTable[["Qx"]],  dplyr::select(yfeat, starts_with("Q_")),
            cTable[["adductx"]],dplyr::select(yfeat, starts_with("adduct_")))
    }
    else if(!is.null(xfeat) & !is.null(yfeat)){
        xfeat <- dplyr::slice(xfeat, rep(seq(1,n()), times = xreps))
        yfeat <- yfeat[yreps,]
        featdata <- cbind.data.frame(
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
    names(featdata) <- paste(rep(fields, each = length(datasets(object))),
                             rep(datasets(object),length(fields)), sep = "_")
    featdata <- cbind.data.frame(rowID = cTable[["rowID"]], featdata)
    return(featdata)
}

#' @title Form combinedTable and featData
#'
#' @description This wraps together the formation of the combinedTable and
#' featdata from the input datasets used to construct the metabCombiner object.
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
    if(any(groupCountX * groupCountY >= 10000))
        stop("irregular group size detected (n > 10000); check m/z values")
    xreps <- rep(groupCountY, times = groupCountX)
    ends <- c(which(!duplicated(yset[["group"]])),nrow(yset)+1)
    yreps <- unlist(lapply(seq(1,nGroups), function(number){
        rep(seq(ends[number], ends[number+1]-1), times = groupCountX[number])
    }))
    cTable <- formCombinedTable(object, xset, yset, xreps, yreps)
    featdata <- form_featdata(object, cTable, xfeat, yfeat, xreps, yreps)
    object <- update_mc(object, combinedTable = cTable, featdata = featdata)
    return(object)
}
