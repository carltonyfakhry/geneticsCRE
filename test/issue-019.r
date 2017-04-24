library(devtools)
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
load_all()

res2 <- GetBestPaths('/data/gcre/segfault-random/randomized_junk_data.txt', nCases = 1537, nControls = 1537, method = 1,
  threshold_percent = 0.05, K = 4, pathLength = 4, iterations = 0, strataF = NA, nthreads = 4)

res2[,c("Paths", "Scores", "Pvalues")]
