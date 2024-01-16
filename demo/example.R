#This is a template script for running metabCombiner

library(metabCombiner)

# read in datasets from file; be sure that stringsAsFactors = FALSE
dataset1 <- read.delim("file/path/to/dataset1.txt", sep = "\t",
                      stringsAsFactors = FALSE)
dataset2 <- read.delim("file/path/to/dataset2.txt", sep = "\t",
                      stringsAsFactors = FALSE)

##column headings in example datasets
names(dataset1) #identity, m/z, rt, adduct, samp1, samp2, samp3, samp4, ...

names(dataset2) #ID, mz, RT, gr1.1, gr1.2, gr1.3, gr2.1, gr2.2, gr2.3,...

############## use metabData to stage each dataset individually ##############
data1 <- metabData(dataset1, mz = "m/z", rt = "rt", id = "identity",
                  adduct = "adduct", samples = "samp", extra = NULL,
                  rtmin = 0.5, rtmax = 29, misspc = 50, measure = "median",
                  zero = TRUE, duplicate = opts.duplicate())

data2 <- metabData(dataset2, mz = "mz", rt = "RT", id = "ID",  adduct = NULL,
                  samples = "gr1", extra = "gr2", rtmin = 0.5, rtmax = 29.5,
                  misspc = 50, measure = "median", zero = TRUE,
                  duplicate = opts.duplicate())

getSamples(data2)   #check that sample names are correct
getExtra(data2)     #check that extra column names are correct
getStats(data2)     #feature statistics

######## Create metabCombiner Object and Group Paired Features by m/z #########

data.combined <- metabCombiner(xdata = data1, ydata = data2, binGap = 0.005,
                              xid = "d1", yid = "d2")

#at any point, this method tracks the progress of the alignment
data.report <- combinedTable(data.combined)

########################### Compute RT Mapping ################################
data.combined <- selectAnchors(data.combined, useID = TRUE, windx = 0.03,
                              windy = 0.03, tolmz = 0.003, tolQ = 0.3)

anchors <- getAnchors(data.combined)   #to view the results of anchor selection

set.seed(100)
data.combined <- fit_gam(data.combined, useID = TRUE, k = seq(12, 20, 2),
                        iterFilter = 2, coef = 2, prop = 0.5, bs = "bs",
                        family = "scat", weights = 1, method = "REML",
                        optimizer = "newton")

##visual of mapping results; see ?plot_fit
plot(data.combined, fit = "gam", main = "Example Fit", xlab = "data1",
     ylab = "data2", pch = 19, lcol = "red", pcol = "black", outlier = "s")

###################### Score Feature Pair Alignments ##########################

# optional function; only run if you have sufficiently representative shared IDs
scores <- evaluateParams(data.combined, A = seq(60, 150, by = 10),
                        B = seq(6, 20), C = seq(0.1, 1 ,0.1), fit = "gam",
                        usePPM = FALSE, minScore = 0.7,
                        penalty = 10, groups = NULL)

data.combined <- calcScores(data.combined, A = 90, B = 15, C = 0.5,
                            fit = "gam", usePPM = FALSE, groups = NULL)


################### Reduce Feature Pair Alignment Report ####################

#option 1 is for routine automated tasks;

#option 1: fully reduced table of 1-1 alignments
data.combined <- reduceTable(data.combined, maxRankX = 2, maxRankY = 2,
                           minScore = 0.5, delta = 0.1)
#this is equivalent to the above reduceTable call
data.combined <- labelRows(data.combined, maxRankX = 2, maxRankY = 2,
                           minScore = 0.5, delta = 0.1, method = "score",
                           resolveConflicts = TRUE, remove = TRUE)


#options 2 and 3 allow for inspection of program determinations
#option 2: fully-detailed row labels, including conflicting feature pair rows
data.combined <- labelRows(data.combined, maxRankX = 2, maxRankY = 2,
                          minScore = 0.5, delta = 0.1, method = "score",
                          resolveConflicts = FALSE, remove = FALSE)

#option 3: resolved conflicting rows, but no row removal
data.combined <- labelRows(data.combined, maxRankX = 2, maxRankY = 2,
                           minScore = 0.5, delta = 0.1, method = "score",
                           resolveConflicts = FALSE, remove = FALSE)

################### Inclusion of Non-Aligned Features ####################

#For including features removed in the previous step or not grouped by m/z
# use the xdata and ydata that were originally used to construct object
data.combined <- updateTables(data.combined, xdata = data1, ydata = data2)

#Aside: updateTables can also be used to make changes to combinedTable report
data.report <- combinedTable(data.combined)
data.report.change <-  dplyr::filter(data.report, mzx > 100 & mzx < 1000)
data.combined <- updateTables(data.combined, combinedTable = data.report.change)

################### Print Combined Feature Table Report ####################

#full table
data.report <- combinedTable(data.combined)
write.table(data.report, file = "combined.dataset.report.csv", sep = ",",
            na = "", row.names = FALSE)

##space-separated by group
write2file(data.report, file = "combined.dataset.report.csv", sep = ",")

######################## Miscellaneous Commands #############################

data.combined     #printed summary of metabCombiner object

getStats(data.combined)  #important object statistics

model <- getModel(data.combined, fit = "gam")  #obtain RT mapping model
getSamples(data.combined, data = "x")   #X dataset sample names
getSamples(data.combined, data = "y")   #Y dataset sample names
getExtra(data.combined, data = "x")     #X dataset extra names
getExtra(data.combined, data = "y")     #Y dataset extra names
getCoefficients(data.combined)   #last used A,B,C weight arguments

