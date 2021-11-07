
#' Resolve Conflicting Alignment Subgroups
#'
#' @description This method resolves conflicting feature pair assignments
#' (labeled as "CONFLICT") to obtain 1-1 feature matches in the
#' \code{combinedTable} results report.
#'
#' @param fields data frame containing the main
#'
#' @param rtOrder logical option to impose RT order for resolving subgroups
#'
#' @details This is called from within \code{\link{labelRows}} (with argument
#' resolveConflicts set to TRUE), \code{\link{reduceTable}}, &
#' \code{\link{metabCombiner}} (using \code{metabCombiner} object inputs). The
#' method determines which combination of unique feature pairs has the highest
#' sum of scores ("resolveScore") within each subgroup. By default, these
#' combinations of feature pairs must have consistency in their retention time
#' order (rtOrder = TRUE). The combination of 1-1 feature pair alignments with
#' the highest resolveScore within the subgroups are annotated as "RESOLVED",
#' with the remaining unannotated rows labeled as "REMOVE" (or removed outright
#' by other package functions). Feature pairs belonging to multiple subgroup
#' (alt > 0) are labeled as REMOVE.
#'
#' @return data.frame of \code{combinedTable} fields, replacing "CONFLICT"
#' labels with "RESOLVED" or "REMOVE", depending on the computations performed.
#'
resolveRows <- function(fields, rtOrder){
    if(is.null(fields[["labels"]]))
        stop("labelRows() has not been called on object")
    if(!any(fields[["labels"]] == "CONFLICT", na.rm = TRUE))
        return(fields)
    fields[["resolveScore"]] <- 0

    confs <- dplyr::arrange(dplyr::filter(fields, .data$subgroup > 0),
                                .data$subgroup, desc(.data$score))
    confs[["labels"]] <- ifelse(confs[["alt"]] > 0, "REMOVE", confs[["labels"]])

    confs[["labels"]] <- .Call("resolveRows", labels = confs[["labels"]],
                                subgroup = confs[["subgroup"]],
                                mzx = confs[["mzx"]], mzy = confs[["mzy"]],
                                rtx = confs[["rtx"]], rty = confs[["rty"]],
                                score = confs[["score"]], rtOrder = rtOrder,
                                resolveScore = confs[["resolveScore"]],
                                PACKAGE = "metabCombiner")

    confs[["labels"]] <- ifelse(confs[["labels"]] == "CONFLICT", "REMOVE",
                                confs[["labels"]])

    confIDs <- match(confs[["rowID"]], fields[["rowID"]])
    fields$labels[confIDs] <- confs[["labels"]]
    fields$resolveScore[confIDs] <- confs[["resolveScore"]]

    return(fields)
}



