
##parameter checks for model fitting functions
check_fit_pars <- function(anchors, fit, useID, iterFilter, ratio, frac,
                            k, spans, iterLoess = 10)
{
    if(nrow(anchors) == 0)
        stop("missing 'anchors' in metabCombiner object; ",
            "selectAnchors() has not been called yet")

    if(nrow(anchors) < 20)
        stop("number of anchors too small (less than 20); ",
            "examine parameters from selectAnchors step.")

    if(!is.logical(useID))
        stop("argument 'useID' must be a logical")

    if((!is.numeric(iterFilter) & !is.integer(iterFilter)) | iterFilter < 0)
        stop("argument 'iterFilter' must be a positive integer")

    if((!is.numeric(ratio) & !is.integer(ratio)) | ratio < 1)
        stop("argument 'ratio' must be a numeric constant greater than 1")

    if((!is.numeric(frac) & !is.integer(frac)) | frac <= 0 | frac > 1)
        stop("argument 'frac' must be numeric between 0 & 1")

    if(fit == "gam"){
        if(!is.numeric(k) & !is.integer(k))
            stop("argument 'k' must be an integer vector")

        if(length(k) == 0)
            stop("argument 'k' must have at least one element")

        if(any(k < 3) | any(k >= nrow(anchors)))
            stop("all values of argument k must be between 3 and n")
    }

    else if(fit == "loess"){
        if(!is.numeric(spans) | any(spans <= 0 | spans >= 1))
            stop("spans argument must be numeric with values between 0 & 1")

        if(length(spans) == 0)
            stop("argument 'spans' must have at least one element")

        if((!is.numeric(iterLoess) & !is.integer(iterLoess)) | iterLoess < 1)
            stop("argument 'iterLoess' must be a >1 integer")
    }

    return(invisible())
}


###scoring parameter checks
check_score_pars <- function(cTable, A, B, C, model, fit, groups,
                            minScore = 0.5, penalty = 5, adduct = 1)
{
    coefs = list(`A` = A, `B` = `B`, `C` = C)

    g = vapply(coefs, function(c) !is.numeric(c) & !is.integer(c), logical(1))
    if(any(g))
        stop("arguments 'A', 'B', 'C' must be numeric constants")

    if(any(A < 0) | any(B < 0) | any(C < 0))
        stop("arguments A, B, C must consist only of non-negative values")

    if(is.null(model))
        stop(paste("object missing model of type ", fit, sep =""))

    if(any(!groups %in% cTable[["group"]]))
        stop("invalid argument 'groups'- one or more groups not found")

    if(!is.numeric(adduct) | adduct < 1)
        stop("argument 'adduct' must be a numeric greater than 1")

    if(!is.numeric(minScore) | minScore < 0 | minScore > 1)
        stop("argument 'minScore' must be a positive constant between 0 & 1")

    if((!is.numeric(penalty) & !is.numeric(penalty)) | penalty < 0)
        stop("argument 'penalty' must be a positive numeric constant")

    return(invisible())
}


##row labeling parameter checks
check_lblrows_pars <- function(maxRankX, maxRankY, minScore, balanced,
                                method, conflict)
{
    if((!is.numeric(maxRankX) & !is.integer(maxRankX)) |
        (!is.numeric(maxRankY) & !is.integer(maxRankY)))
        stop("arguments maxRankX & maxRankY must be integers")

    if(maxRankX < 1 | maxRankY < 1)
        stop("arguments maxRankX & maxRankY must be equal to 1 or greater")

    if(!is.numeric(minScore) | minScore >= 1 | minScore <= 0)
        stop("argument 'minScore' must be a numeric value between 0 and 1")

    if(!is.logical(balanced))
        stop("expected a logical for argument 'balanced'")

    if(method == "score"){
        if(!is.numeric(conflict) | conflict > 1 | conflict < 0)
            stop("argument 'conflict' must be a numeric value between 0 & 1")

        if(length(conflict) > 1)
            stop("constant value expected for conflict")
    }

    else if(method == "mzrt"){
        if(length(conflict) != 4)
            stop("with method = mzrt, length 4 vector expected for conflict)")

        if(any(conflict < 0))
            stop("values in 'conflict' argument must be non-negative")
    }

    return(invisible())
}



