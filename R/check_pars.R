
#parameter checks for metabCombiner() constructor step
check_combine_pars <- function(binGap, means, xid, yid)
{
    if(is.null(binGap) | length(binGap) == 0)
        stop("'binGap' argument must be a defined constant")
    if(binGap <= 0 | !is.numeric(binGap))
        stop("'binGap' argument must be a positive numeric constant")
    if(binGap >= 0.1)
        stop("'binGap' argument is too large; recommended range: 0.005-0.01")
    if(binGap > 0.05)
        warning("large 'binGap' argument value; recommended range: 0.005-0.01")
    if(!is.logical(unlist(means)))
        stop("arguments means must contain only logicals")
}


#parameter checks for anchor selection
check_anchors_pars <- function(useID, tolmz, tolQ, tolrtq, windx, windy)
{
    pars <- list(tolmz, tolQ, tolrtq, windx, windy)
    if(any(vapply(pars, function(p) !is.numeric(p), logical(1))))
        stop("arguments tolmz, tolQ, tolrtq, windx & windy must be numeric")

    if(any(vapply(pars, function(p) p <= 0, logical(1))))
        stop("arguments tolmz, tolQ, tolrtq, windx & windy must be positive")

    if(!is.logical(useID)) stop("expected logical value for argument 'useID'")
}


##parameter checks for model fitting functions
check_fit_pars <- function(anchors, fit, useID, iterFilter, coef, prop, k, spans)
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

    if((!is.numeric(coef) & !is.integer(coef)) | coef < 1)
        stop("argument 'coef' must be a numeric constant greater than 1")

    if((!is.numeric(prop) & !is.integer(prop)) | prop <= 0 | prop > 1)
        stop("argument 'prop' must be numeric between 0 & 1")

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
    }
}


###scoring parameter checks
check_score_pars <- function(cTable, A, B, C, model, fit, groups,
                            minScore = 0.5, penalty = 5, adduct = 1)
{
    coefs <- list(`A` = A, `B` = `B`, `C` = C)

    g = vapply(coefs, function(c) !is.numeric(c) & !is.integer(c), logical(1))
    if(any(g))
        stop("arguments 'A', 'B', 'C' must be numeric constants")

    if(any(A < 0) | any(B < 0) | any(C < 0))
        stop("arguments A, B, C must consist only of non-negative values")

    if(is.null(model))
        stop(paste("object missing model of type ", fit, sep =""))

    if(any(!groups %in% cTable[["group"]]))
        stop("invalid argument 'groups'- one or more groups not found")

    if(any(groups <= 0))
        stop("'groups' argument must be strictly positive")

    if(!is.numeric(adduct) | adduct < 1)
        stop("argument 'adduct' must be a numeric greater than 1")

    if(!is.numeric(minScore) | minScore < 0 | minScore > 1)
        stop("argument 'minScore' must be a positive constant between 0 & 1")

    if((!is.numeric(penalty) & !is.numeric(penalty)) | penalty < 0)
        stop("argument 'penalty' must be a positive numeric constant")
}


##row labeling parameter checks
check_lblrows_pars <- function(maxRankX, maxRankY, minScore, maxRTerr, balanced,
                                method, delta, resolveConflicts, rtOrder)
{
    if((!is.numeric(maxRankX) & !is.integer(maxRankX)) |
        (!is.numeric(maxRankY) & !is.integer(maxRankY)))
        stop("arguments maxRankX & maxRankY must be integers")

    if(maxRankX < 1 | maxRankY < 1)
        stop("arguments maxRankX & maxRankY must be equal to 1 or greater")

    if(!is.numeric(minScore) | minScore >= 1 | minScore <= 0)
        stop("argument 'minScore' must be a numeric value between 0 and 1")

    if(!is.numeric(maxRTerr) | maxRTerr <= 0)
        stop("argument 'maxRTerr' must be a numeric value greater than 0")

    if(!is.logical(balanced))
        stop("expected a logical for argument 'balanced'")

    if(!is.logical(resolveConflicts))
        stop("expected a logical for argument 'resolveConflicts'")

    if(!is.logical(rtOrder))
        stop("expected a logical for argument 'rtOrder'")

    if(method == "score"){
        if(!is.numeric(delta) | delta > 1 | delta < 0)
            stop("argument 'delta' must be a numeric value between 0 & 1")
        if(length(delta) > 1)
            stop("constant value expected for argument 'delta'")
    }

    else if(method == "mzrt"){
        if(length(delta) != 4)
            stop("with method = 'mzrt', length 4 vector expected for delta)")

        if(!is.numeric(delta) | any(delta < 0))
            stop("values in 'delta' argument must be non-negative")
    }
}



