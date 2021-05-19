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
#' @param pcol color of the normal points (ordered RT pair) in the plot
#'
#' @param pch plot character type; see ?graphics::par for details
#'
#' @param lcol color of the fitted line in the plot
#'
#' @param lwd line width of the curve fit between anchor points
#'
#' @param outlier display option for outliers. If "show" or "s", treats
#' outlier points like normal anchors; if "remove" or "r", removes outlier
#' points from the plot; if "highlight" or "h", displays outliers with a
#' different color and associated legend.
#'
#' @param ocol color of the outlier points; outlier argument must be set to
#' "highlight" or "h"
#'
#' @param legend length-2 character vector indicating point labels in the legend
#' if outlier argument set to "highlight" or "h"
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
#' p.comb = fit_gam(p.comb, k = 20, iterFilter = 1, family = "gaussian")
#'
#' ##plot of GAM fit
#' plot(p.comb, main = "Example GAM Fit Plot", xlab = "X Dataset RTs",
#'      ylab = "Y Dataset RTs", pcol = "red", lcol = "blue", lwd = 5,
#'      fit = "gam", outliers = "remove")
#'
#' grid(lwd =  2, lty = 3 ) #adding gridlines
#'
#' @export
plot_fit <- function(object, fit = c("gam","loess"), pcol = "black",
                    lcol = "red", lwd = 3, pch = 19, outlier = "show",
                    ocol = "springgreen4", legend = c("anchor", "outlier"),
                    ...)
{
    fit <- match.arg(fit)
    model <- getModel(object, fit = fit)

    if(is.null(model))
        stop(paste("object missing model of type", fit))

    if(length(legend) != 2)
        stop("legend must be a length 2 character vector")

    if(fit == "loess")
        data <- data.frame(rtx = model[["x"]], rty = model[["y"]],
                        pred = model[["fitted"]], weights = model[["weights"]])

    else if (fit == "gam")
        data <- data.frame(rtx = model$model$rtx, rty = model$model$rty,
                            preds = model[["fitted.values"]],
                            weights = model[["prior.weights"]])

    data <- dplyr::arrange(data, .data$rtx)

    if(outlier %in% c("remove", "r", "highlight", "h")){
        data2 <- dplyr::filter(data, .data$weights == 0)
        data <- dplyr::filter(data, .data$weights > 0)
    }

    rtx <- data[["rtx"]];  rty <- data[["rty"]]
    graphics::plot(rtx, rty, type = "p", col = pcol, pch = pch, ...)
    graphics::lines(x = data[["rtx"]], y = data[["preds"]],
                    col = lcol, lwd = lwd,...)

    if(outlier %in% c("highlight", "h")){
        graphics::points(data2[["rtx"]],data2[["rty"]], col = ocol,
                        pch = pch,...)
        graphics::legend("topleft", legend = legend, col = c(pcol, ocol),
                        pch = pch)
    }
}
