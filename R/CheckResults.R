#' This function checks the results returned by GetBestPaths
#' @export
checkBestPaths <- function(dataset, BestPaths, pathLength, nCases, nControls, method){

  # Precompute values table
  print("Computing Values Table...")
  ValueTable <- geneticsCRE:::getValuesTable(nCases, nControls)

  # Preprocessing Phenotype dataset, all genes which occur in more
  # than threshold_percent of patients will be removed from the
  # dataset, also genes which do not occur in any patients will be removed
  print("Preprocessing Phenotype dataset...")
  threshold_percent <- 1
  dataset <- geneticsCRE:::PreprocessTable(dataset, threshold_percent, nCases, nControls)
  genes <- dataset$genes
  data <- dataset$data
  genes_data <- data.frame(genes = genes, data, stringsAsFactors = F)
  patients <- dataset$patients
  rm(dataset)

  passed <- TRUE

  # Handle paths of length 2 and greater
  for(path_len in 1:pathLength){

    bestpaths <- BestPaths[BestPaths$Lengths == path_len,]

    for(i in 1:nrow(bestpaths)){

      # Split a path which is a string of the form "Gene1 (-) -> Gene2 (+)" as an example
      split_path <- unlist(strsplit(bestpaths$SignedPaths[i], split = " -> "))

      # Extract the signs from the split path
      sign <- sapply(split_path, function(x) unlist(strsplit(x, split = " "))[2])
      sign <- ifelse(sign == "(+)", 1, -1)

      # Get the genes in the path
      genes <- sapply(split_path, function(x) unlist(strsplit(x, split = " "))[1])

      # Get the data for each gene in the path
      path_data_pos <- as.matrix(genes_data[match(genes, genes_data$genes),2:ncol(genes_data)],drop = FALSE)
      path_data_neg <- matrix(0, nrow = nrow(path_data_pos), ncol = ncol(path_data_pos))
      if(method == 2){
        path_data_neg[which(sign == -1),] = path_data_pos[which(sign == -1),]
        path_data_pos[which(sign == -1),] = 0
      }

      flipped <- ifelse(sign[1] == -1, T, F)
      names(flipped) <- NULL

      subpath_pos1 <- colSums(path_data_pos)
      subpath_pos1[subpath_pos1 != 0] <- 1
      subpath_neg1 <- colSums(path_data_neg)
      subpath_neg1[subpath_neg1 != 0] <- 1

      inds_pos1 <- which(subpath_pos1 != 0)
      inds_neg1 <- which(subpath_neg1 != 0)
      # intersect_inds <- intersect(inds_pos1, inds_neg1)
      # inds_pos1 <- setdiff(inds_pos1, intersect_inds)
      # inds_neg1 <- setdiff(inds_neg1, intersect_inds)

      # compute the score of the path
      cases_pos <- length(which(inds_pos1 <= nCases))
      controls_pos <- length(inds_pos1) - cases_pos
      cases_neg <- length(which(inds_neg1 > nCases))
      controls_neg <- length(inds_neg1) - cases_neg
      cases <- cases_pos + cases_neg
      controls <- controls_pos + controls_neg
      # score <- ifelse(flipped & method == 1, ValueTable[(controls+1),(cases+1)], ValueTable[(cases+1),(controls+1)])
      if(method == 1){
        score <- ValueTable[(cases+1),(controls+1)]
      }else if(method == 2){
        score <- ValueTable[(cases_pos+1),(controls_pos+1)] + ValueTable[(cases_neg+1),(controls_neg+1)]
      }

      # Check for score equality
      if(score != bestpaths$Scores[i]){
        passed <- FALSE
        print(sprintf("bad score for: '%s'", bestpaths[i,1]))
        print(sprintf("cases: %d, ctrls: %d, check score: %f", cases, controls, score))
        print(bestpaths[i,])
      }

    }

  }

  return(passed)

}
