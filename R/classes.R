##class descriptions go in here
#

############################# Classes #####################

#'
#' @exportClass
setClass("metabData", slots = c(data = "data.frame",
                                samples = "character",
                                extra = "character",
                                stats = "list"),
                      prototype = prototype(
                                 data = data.frame(),
                                 samples = character(),
                                 extra = character(),
                                 stats = list()
                      )
)

#'
#' @exportClass
setClass("metabCombiner", slots = c(xdata = "metabData",
                                    ydata = "metabData",
                                    combinerTable = "data.frame",
                                    binGap = "numeric",
                                    anchors = "data.frame",
                                    model = "list",
                                    coefficients = "list",
                                    stats = "list"
                          ),
                          prototype = prototype(
                                   xdata = new("metabData"),
                                   ydata = new("metabData"),
                                   combinerTable = data.frame(),
                                   binGap = numeric(),
                                   anchors = data.frame(),
                                   model = list(),
                                   coefficients = list(),
                                   stats = list()
                          )
)




