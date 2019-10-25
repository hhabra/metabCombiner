#' @title Plot Combiner Fits
#'
#' @description
#' This is the plotting method for Combiner objects. Displays anchoring ordered
#' pairs and a curve fit computed using \code{fit_gam} or \code{fit_loess}.
#'
#' @param object metabCombiner object
#'
#' @param fit Choice of model (either "gam" or "loess").
#'
#' @param pcol color of the points (ordered pairs) in the plot.
#'
#' @param lcol color of the fitted line in the plot
#'
#' @param lwd line width of the curve fit between anchor points.
#'
#' @param ... Other variables passed into graphics::plot
#'
#' @export
plot_Combiner <- function(object, fit = c("gam","loess"), pcol, lcol, lwd,...)
{
    fit = match.arg(fit)
    model = getModel(object, fit = fit)

    if(is.null(model))
        stop(paste("object missing model of type", fit))

    if(fit == "loess")
        data = data.frame(rtx = model[["x"]],
                          rty = model[["y"]],
                          pred = model[["fitted"]],
                          weights = model[["weights"]])

    else if (fit == "gam")
        data = data.frame(rtx = model$model$rtx,
                          rty = model$model$rty,
                          pred = model[["fitted.values"]],
                          weights = model[["prior.weights"]])

    data = dplyr::arrange(data, rtx)

    if(missing(pcol))
        pcol = "black"
    if(missing(lcol))
        lcol = "red"
    if(missing(lwd))
        lwd = 4

    graphics::plot(data[["rtx"]], data[["rty"]], type = "p", col = pcol,...)
    graphics::lines(x = data[["rtx"]], y = data[["pred"]],
                    col = lcol, lwd = lwd,...)
}

