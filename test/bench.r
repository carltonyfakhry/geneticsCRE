library(devtools)
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
load_all()

# res2 <- GetBestPaths('/data/gcre/testdata', nCases = 1537, nControls = 1537, method = 'method2-old', threshold_percent = 0.05, K = 2, pathLength = 5, iterations = 4, strataF = NA, nthreads = 2)
# res2 <- GetBestPaths('/data/gcre/testdata2', nCases = 1000, nControls = 1000, method = 'method2', threshold_percent = 0.05, K = 8, pathLength = 3, iterations = 200, strataF = NA, nthreads = 4)
# res2 <- GetBestPaths('/data/gcre/data_null', nCases = 100, nControls = 100, method = 'directed', threshold_percent = 0.05, K = 4, pathLength = 3, iterations = 10, strataF = NA, nthreads = 4)
# res2 <- GetBestPaths('/data/sim_lof_500', nCases = 500, nControls = 500, method = 'enrich', threshold_percent = 1, K = 8, pathLength = 3, iterations = 10, strataF = NA, nthreads = 4)

res2 <- GetBestPaths('/data/gcre/testdata', nCases = 1537, nControls = 1537, method = 1, threshold_percent = 0.05, K = 4, pathLength = 4, iterations = 1000, strataF = NA, nthreads = 4)
res2[,c("Paths", "Scores", "Pvalues")]

res2 <- GetBestPaths('/data/gcre/testdata', nCases = 1537, nControls = 1537, method = 2, threshold_percent = 0.05, K = 4, pathLength = 4, iterations = 1000, strataF = NA, nthreads = 4)
res2[,c("Paths", "Scores", "Pvalues")]

res2 <- GetBestPaths('/data/gcre/sim_lof_500', nCases = 500, nControls = 500, method = 'enrich', threshold_percent = 1, K = 4, pathLength = 5, iterations = 100, strataF = NA, nthreads = 4)
res2[,c("Paths", "Scores", "Pvalues")]
