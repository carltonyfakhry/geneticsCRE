library(devtools)
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
load_all()

data <- '/data/gcre/testdata'
len <- 1
m <- 2

res2 <- GetBestPaths(dataset = data, nCases = 1537, nControls = 1537, method = m,
threshold_percent = 1, K = 1000, pathLength = len, iterations = 10, strataF = NA, nthreads = 4)

checkBestPaths(dataset = data, BestPaths = res2, pathLength = len, nCases = 1537, nControls = 1537, method = m)

# res2[,c("Paths", "Scores", "Pvalues")]
