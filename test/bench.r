library(devtools)
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
load_all()

res2 <- GetBestPaths('/data/gcre/testdata', nCases = 1537, nControls = 1537, method = 1, threshold_percent = 1, K = 12, pathLength = 4, iterations = 10, strataF = NA, nthreads = 4)

res2[,c("Scores", "Pvalues", "Cases", "Controls", "debug")]
