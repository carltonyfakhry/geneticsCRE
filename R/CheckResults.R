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

      subpath_pos2[intersect(inds_pos1, inds_pos2)] <- 0
      subpath_neg2[intersect(inds_neg1, inds_neg2)] <- 0
      inds_pos2 <- setdiff(inds_pos2, inds_pos1)
      inds_neg2 <- setdiff(inds_neg2, inds_neg1)

      # compute the score of the path
      cases1 <- sum(subpath_pos1[1:nCases]) + sum(subpath_neg1[nCases:(nCases+nControls)])
      controls1 <- sum(subpath_pos1[nCases:(nCases+nControls)]) + sum(subpath_neg1[1:nCases])
      cases2 <- sum(subpath_pos2[1:nCases]) + sum(subpath_neg2[nCases:(nCases+nControls)])
      controls2 <- sum(subpath_pos2[nCases:(nCases+nControls)]) + sum(subpath_neg2[1:nCases])
      score <- ifelse(flipped & method == 1, ValueTable[(controls1+controls2+1),(cases1+cases2+1)], ValueTable[(cases1+cases2+1),(controls1+controls2+1)])

      # Check for score equality
      if(score != bestpaths$Scores[i]){
        stop(paste("The following path:", bestpaths[i,1], "has the wrong score!", sep = " "))
      }

    }

  }

  return(TRUE)

}
