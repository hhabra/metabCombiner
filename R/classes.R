##class descriptions go in here
#

############################# Classes ##################### 

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

setClass("metabCombiner", slots = c(xdata = "metabData",
                                    ydata = "metabData",
                                    binGap = "numeric",
                                    anchors = "data.frame",
                                    model = "list",
                                    scores = "list",
                                    combinerTable = "data.frame"
                          ),
                          prototype = prototype(
                                   xdata = new("metabData"),
                                   ydata = new("metabData"),
                                   binGap = numeric(),
                                   anchors = data.frame(),
                                   model = list(),
                                   scores = list(),
                                   combinerTable = data.frame()
                          )
)




