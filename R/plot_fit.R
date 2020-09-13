#' @title Plot metabCombiner Fits
#'
#' @description
#' This is a plotting method for metabCombiner objects. It displays
#' ordered pairs and a curve fit computed using \code{fit_gam} or
#' \code{fit_loess}, using base R graphics.
#'
#' @param object metabCombiner object
#'
#' @param fit choice of model (either "gam" or "loess").
#'
#' @param pcol color of the points (ordered pairs) in the plot.
#'
#' @param lcol color of the fitted line in the plot
#'
#' @param lwd line width of the curve fit between anchor points
#'
#' @param remove.outliers logical, option to exclude outlier anchors in plot
#'
#' @param ... Other variables passed into graphics::plot
#'
#' @return no values returned
#'
#' @examples
#' data(plasma30)
#' data(plasma20)
#'
#' p30 <- metabData(plasma30, samples = "CHEAR")
#' p20 <- metabData(plasma20, samples = "Red", rtmax = 17.25)
#' p.comb = metabCombiner(xdata = p30, ydata = p20, binGap = 0.0075)
#' p.comb = selectAnchors(p.comb, tolmz = 0.003, tolQ = 0.3, windy = 0.02)
#' p.comb = fit_gam(p.comb, k = 20, iterFilter = 1)
#'
#' ##plot of GAM fit
#' plot(p.comb, main = "Example GAM Fit Plot", xlab = "X Dataset RTs",
#'      ylab = "Y Dataset RTs", pcol = "red", lcol = "blue", lwd = 5,
#'      fit = "gam", pch = 19, remove.outliers = TRUE)
#'
#' grid(lwd =  2, lty = 3 ) #adding gridlines
#'
#' @export
plot_fit <- function(object, fit = c("gam","loess"), pcol, lcol, lwd,
                    remove.outliers = FALSE, ...)
{
    fit = match.arg(fit)
    model = getModel(object, fit = fit)

    if(is.null(model))
        stop(paste("object missing model of type", fit))

    if(fit == "loess")
        data = data.frame(rtx = model[["x"]], rty = model[["y"]],
                            pred = model[["fitted"]],
                            weights = model[["weights"]])

    else if (fit == "gam")
        data = data.frame(rtx = model$model$rtx, rty = model$model$rty,
                            preds = model[["fitted.values"]],
                            weights = model[["prior.weights"]])

    data = dplyr::arrange(data, .data$rtx)

    if(remove.outliers)
        data = dplyr::filter(data, .data$weights > 0)

    if(missing(pcol))
        pcol = "black"
    if(missing(lcol))
        lcol = "red"
    if(missing(lwd))
        lwd = 4

    rtx = data[["rtx"]]
    rty = data[["rty"]]

    graphics::plot(rtx, rty, type = "p", col = pcol,...)
    graphics::lines(x = data[["rtx"]], y = data[["preds"]],
                    col = lcol, lwd = lwd,...)
}
