library(devtools)
load_all()

# res2 <- GetBestPaths('/data/gcre/testdata', nCases = 1537, nControls = 1537, method = 'method2-old', threshold_percent = 0.05, K = 2, pathLength = 5, iterations = 4, strataF = NA, nthreads = 2)
# res2 <- GetBestPaths('/data/gcre/testdata', nCases = 1537, nControls = 1537, method = 'method2-old', threshold_percent = 0.05, K = 2, pathLength = 4, iterations = 80, strataF = NA, nthreads = 1)
res2 <- GetBestPaths('/data/gcre/testdata', nCases = 1537, nControls = 1537, method = 'method2', threshold_percent = 0.05, K = 2, pathLength = 1, iterations = 4, strataF = NA, nthreads = 1)
res2