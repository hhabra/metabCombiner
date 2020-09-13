#This is a template script for running metabCombiner

library(metabCombiner)

# read in datasets from file; be sure that stringsAsFactors = FALSE
dataset1 = read.delim("file_path_to_dataset1.txt", sep = "\t",
                      stringsAsFactors = FALSE)

dataset2 = read.delim("file_path_to_dataset1.txt", sep = "\t",
                      stringsAsFactors = FALSE)


##column headings in example datasets
names(dataset1) #identity, m/z, rt, adduct, samp1, samp2, samp3, samp4, ...

names(dataset2) #ID, mz, RT, gr1.1, gr1.2, gr1.3, gr2.1, gr2.2, gr2.3,...


############## use metabData to stage each dataset individually ##############
data1 = metabData(dataset1, mz = "m/z", rt = "rt", id = "identity",
                  adduct = "adduct", samples = "samp", extra = NULL,
                  rtmin = 0.5,rtmax = , misspc = 50, measure = "median",
                  zero = TRUE, duplicate = c(0.0025,0.05))

data2 = metabData(dataset2, mz = "mz", rt = "RT", id = "ID",  adduct = NULL,
                  samples = "gr1", extra = "gr2", rtmin = 0.5, rtmax = 29.5,
                  misspc = 50, measure = "median", zero = TRUE,
                  duplicate = c(0.0025,0.05))

getSamples(data1)   #check that sample names are correct
getExtra(data1)     #check that extra column names are correct
getStats(data1)     #feature statistics

######## Create metabCombiner Object and Group Paired Features by m/z #########

data.combined = metabCombiner(xdata = data1, ydata = data2, binGap = 0.005)

data.report = combinedTable(data.combined) #template of the final report table

########################### Compute RT Mapping ################################

data.combined = selectAnchors(data.combined, useID = TRUE, windx = 0.03,
                              windy = 0.03, tolmz = 0.003, tolQ = 0.3)

anchors = getAnchors(data.combined)   #to view the results of anchor selection

set.seed(100)
data.combined = fit_gam(data.combined, useID = TRUE, k = seq(12, 20, 2),
                        iterFilter = 2,ratio = 2, frac = 0.5, bs = "bs",
                        family = "scat", weights = 1, method = "REML",
                        optimizer = "newton")

##visual of mapping results
plot(data.combined, fit = "gam", main = "Example Fit", xlab = "data1",
     ylab = "data2", pch = 19, lcol = "red", pcol = "black")

###################### Score Feature Pair Alignments ##########################

# optional function; only run if you have sufficiently representative shared IDs
scores = evaluateParams(data.combined, A = seq(60, 150, by = 10),
                        B = seq(6, 20), C = seq(0.1, 1 ,0.1), fit = "gam",
                        PPM = FALSE, useAdduct = FALSE, minScore = 0.7,
                        penalty = 10, groups = NULL)


data.combined = calcScores(data.combined, A = 100, B = 15, C = 0.5,
                           fit = "gam", usePPM = FALSE, groups = NULL)


################### Reduce Feature Pair Alignment Report ####################

data.report = combinedTable(data.combined)

#new columns added to report table: program-determined labels
data.report = labelRows(data.report, maxRankX = 3,maxRankY = 3, minScore = 0.5,
                        conflict = 0.1, method = "score", balanced = TRUE,
                        remove = FALSE)

################### Print Space-Separated Groups Report ####################

write2file(data.report, file = "combined.dataset.report.csv", sep = ",")


######################## Miscellaneous Commands #############################

data.combined     #printed summary of metabCombiner object

getStats(data.combined)  #important object stats

model = getModel(data.combined, fit = "gam")  #obtain RT mapping model

getSamples(data.combined, data = "x")   #X dataset sample names
getSamples(data.combined, data = "y")   #Y dataset sample names
getExtra(data.combined, data = "x")     #X dataset extra names
getExtra(data.combined, data = "y")     #Y dataset extra names

nonMatchedX = nonmatched(data.combined, data = "x") #non-matched X features

getCoefficients(data.combined)   #last used A,B,C weight arguments

