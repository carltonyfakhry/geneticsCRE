library(devtools)
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
load_all()

# data <- '/data/gcre/issue-020/testdata'
data <- '/data/gcre/testdata2'
len <- 1
m <- 1
cases <- 1000

res2 <- GetBestPaths(dataset = data, nCases = cases, nControls = cases, method = m,
threshold_percent = 1, K = 1000, pathLength = len, iterations = 10, strataF = NA, nthreads = 4)

checkBestPaths(dataset = data, BestPaths = res2, pathLength = len, nCases = cases, nControls = cases, method = m)

# res2[,c("Paths", "Scores", "Pvalues")]
