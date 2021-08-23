
#' Resolve Conflicting Alignment Subgroups
#'
#' @param fields data frame containing the main
#'
#' @param rtOrder logical option to impose RT order for resolving subgroups
#'
#' @noRd
resolveRows <- function(fields, rtOrder){
    if(!any(fields[["labels"]] == "CONFLICT"))
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



