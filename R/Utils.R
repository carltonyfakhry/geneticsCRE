# Randomize columns of data
getRandomizedDatCols <- function(ncol_data, strataF, nCases, nControls, stratagroups, stratanumbers){

  randomized_data_indices <- NA
  if(length(strataF) == 1 & is.na(strataF)){
    randomized_data_indices <- sample(1:ncol_data, replace = F)
  }else{
    randomized_data_indices <- rep(0, ncol_data)
    for(i in stratagroups){
      stratainds <- which(stratanumbers == i)
      samplestratainds <- sample(stratainds, size = length(stratainds), replace = F)
      randomized_data_indices[stratainds] <- samplestratainds
    }
  }
  return(randomized_data_indices)

}



# Get the matrix whose rows are randomized columns of the data
getRandIndicesMat <- function(ncol_data, iterations, strataF, nCases, nControls, stratagroups, stratanumbers){

  RandIndicesMat <- matrix(0, nrow = iterations, ncol = ncol_data)

  if(length(strataF) == 1 & is.na(strataF)){

    for(i in 1:iterations){
      randomized_cols <- geneticsCRE:::getRandomizedDatCols(ncol_data, strataF, nCases, nControls, NA, NA)
      randomized_cols <- randomized_cols - 1
      RandIndicesMat[i,] <- randomized_cols
    }

  }else{

    for(i in 1:iterations){
      randomized_cols <- geneticsCRE:::getRandomizedDatCols(ncol_data, strataF, nCases, nControls, stratagroups, stratanumbers)
      randomized_cols <- randomized_cols - 1
      RandIndicesMat[i,] <- randomized_cols
    }

  }

  return(RandIndicesMat)

}


# Convert uids to protein symbols in paths
Uid2Symbol <- function(Ents, UserKpaths){

  splitted <- strsplit(as.character(UserKpaths), split = "->")
  symbpaths <- rep("", length(splitted))

  for(i in 1:length(splitted)){

    uids <- as.integer(splitted[[i]])
    inds <- match(uids, Ents$uid)
    symbs <- Ents$symbol[inds]
    symbpaths[i] <- paste(symbs, collapse = " -> ")

  }

  return(symbpaths)

}



# Convert signed uids to signed protein symbols in paths
SignUid2SignSymbol <- function(Ents, UserKsignpaths){

  splitted <- strsplit(as.character(UserKsignpaths), split = " -> ")
  symbpaths <- rep("", length(splitted))

  for(i in 1:length(splitted)){

    uids <- c()
    signs <- c()
    for(j in splitted[[i]]){ # j has the form: uid (sign)
      words <- strsplit(j, split = " ")[[1]]
      uids <- c(uids, words[1])
      signs <- c(signs, words[2])
    }

    uids <- as.integer(uids)
    inds <- match(uids, Ents$uid)
    symbs <- Ents$symbol[inds]
    symbpaths[i] <- paste(paste(symbs, signs, sep = " "), collapse = " -> ")

  }

  return(symbpaths)

}


# Prepare data to be in the order in which it will be joined to its corresponding path.
# trgs_toB_joined is a vector of uids whose data will be joined respectively
# to some path.
matchData <- function(genes_data, trgs_toB_joined){

  genes_data <- genes_data[match(trgs_toB_joined, genes_data$uid),]
  data <- as.matrix(genes_data[,3:ncol(genes_data)])
  return(data)

}


matchDataIndices <- function(genes_data, trgs_toB_joined){

  matched_indices <- match(trgs_toB_joined, genes_data$uid)
  return(matched_indices)

}



# Get a list of empty paths
getListEmptyPaths <- function(len_counts){

  lst <- list()
  lst[["counts_pos"]] <- rep(0, len_counts)
  lst[["counts_neg"]] <- rep(0, len_counts)
  lst[["counts_conflict"]] <- rep(0, len_counts)
  lst[["cols_pos"]] <- vector(mode="numeric", length=0)
  lst[["cols_neg"]] <- vector(mode="numeric", length=0)
  lst[["cols_conflict"]] <- vector(mode="numeric", length=0)
  lst[["locations_pos"]] <- rep(0, len_counts)
  lst[["locations_neg"]] <- rep(0, len_counts)
  lst[["locations_conflict"]] <- rep(0, len_counts)
  return(lst)

}



# # This function computes the values table
# getValuesTable <- function(nControls, nCases){
#
#   probTable <- matrix(NA, nCases + 1, nCases + nControls + 1)
#   valuesTable <- matrix(NA, nCases + 1, nControls + 1)
#   n <- nCases + nControls
#   for (i in 0:n) {
#     i_vector = (max(0, (i - nControls)):min(i, nCases))
#     probTable[i_vector + 1, i + 1] = dhyper(i_vector, nCases, nControls, i)
#   }
#
#   ## Now compute the p-values
#   for (i in 0:n) {
#     i_max = max(0, (i - nControls))
#     i_min = min(i, nCases)
#     i_vector =  i_max:i_min
#     valuesTable[cbind(i_vector + 1, i - i_vector + 1)] = -log(rev(cumsum(probTable[(i_min:i_max) + 1, i + 1])))
#   }
#   valuesTable[is.infinite(valuesTable)] = max(valuesTable[is.finite(valuesTable)]) + 1;
#   return(valuesTable)
#
# }



# This function computes the values table
getValuesTable <- function(nControls, nCases){

  probTable <- matrix(NA, nCases + 1, nCases + nControls + 1)
  valuesTable <- matrix(NA, nCases + 1, nControls + 1)
  n <- nCases + nControls
  for (i in 0:n) {
    i_vector = (max(0, (i - nControls)):min(i, nCases))
    probTable[i_vector + 1, i + 1] = dhyper(i_vector, nCases, nControls, i)
  }

  ## Now compute the p-values
  for (i in 0:n) {
    i_max = max(0, (i - nControls))
    i_min = min(i, nCases)
    i_vector =  i_max:i_min
    valuesTable[cbind(i_vector + 1, i - i_vector + 1)] = -log(rev(cumsum(probTable[(i_min:i_max) + 1, i + 1])))
  }
  valuesTable[is.infinite(valuesTable)] = max(valuesTable[is.finite(valuesTable)]) + 1;
  return(valuesTable)

}



# Preprocess phenotype data
PreprocessTable <- function(dataset, threshold_percent, nCases, nControls){

  df <- read.table(dataset, sep = "\t", header = T, stringsAsFactors = F)

  if(names(df)[1] != "symbols"){
    stop("First column must be named symbols and must contain the gene symbols!")
  }

  df <- df[!is.na(df$symbols),]
  df <- df[!duplicated(df$symbols),]
  df[df == 2] <- 1

  if((nCases + nControls) != (ncol(df)-1)){
    stop("The number of patients in the dataset must be equal to nCases + nControls!")
  }

  if(any(df[,2:ncol(df)] != 0 & df[,2:ncol(df)] != 1)){
    stop("The patient columns must consist of only 0,1 or 2 entries!")
  }

  # Keep genes with frequencies less than threshold_percent
  freqs <- rowSums(df[,2:ncol(df)], na.rm = T)
  targetFreq <- threshold_percent * ncol(df)
  inds <- which(freqs <= targetFreq)
  df <- df[inds,]

  genes <- as.vector(df[,1])
  patients <- colnames(df[,2:ncol(df)])
  data <- as.matrix(df[,2:ncol(df)])
  rm(df)
  results <- list(data = data, genes = genes, patients = patients)
  return(results)

}



# Check the input parameters
check_input <- function(nCases, nControls, method, threshold_percent, K, pathLength, iterations, nthreads){

  if(!is.numeric(nCases) | length(nCases) != 1 | floor(nCases) != nCases | nCases < 1){
    stop("nCases must be an integer greater than 0!")
  }

  if(!is.numeric(nControls) | length(nControls) != 1 | floor(nControls) != nControls | nControls < 1){
    stop("nControls must be an integer greater than 0!")
  }

  if(!is.numeric(K) | length(K) != 1 | floor(K) != K | K < 1){
    stop("K must be an integer greater than 0!")
  }

  # methods can have additional tags to specify an implementation (eg 'method2-lomem') --dmitri
  if(!grepl('^method[12]', method)) {
    stop("method must be one of: 'method1', 'method2'!")
  }

  if(!is.numeric(threshold_percent) | length(threshold_percent) != 1 | threshold_percent > 1 | threshold_percent < 0){
    stop("threshold_percent must be a real number greater than or equal to 0 and less than or equal to 1!")
  }

  if(!is.numeric(pathLength) | length(pathLength) != 1 | floor(pathLength) != pathLength | pathLength > 5 | pathLength < 1){
    stop("pathLength must be an integer greater than 0 and less than 6!")
  }

  if(!is.numeric(iterations) | length(iterations) != 1 | floor(iterations) != iterations | iterations < 1){
    stop("iterations must be an integer greater than 0!")
  }

  if((nCases + nControls) > 65536) stop("Package cannot process data set with more than 65536 columns!")

  if(length(nthreads) != 1 | (!is.na(nthreads) & (!is.numeric(nthreads) | nthreads != floor(nthreads)))) stop("nthreads must be set to a postive integer, otherwise should be set to NA!")

  # if(!dir.exists(path)) stop("Please provide a valid path that the package can use to write and read data!")

}



# For columns <= nCases will be set to 1 if they are permuted and remain <= nCases.
# For columns > nControls will be set to 1 if they are permuted and remain > nCases.
getCaseORControl <- function(RandIndicesMat, nCases, nControls){

  CaseORControl <- matrix(0, nrow = nrow(RandIndicesMat), ncol = ncol(RandIndicesMat))

  for(i in 1:nrow(CaseORControl)){
    for(j in 1:nCases){
      CaseORControl[i,j] = ifelse(RandIndicesMat[i,j] < nCases, 1, 0)
    }
    for(j in (nCases+1):(nCases+nControls)){
      CaseORControl[i,j] = ifelse(RandIndicesMat[i,j] >= nCases, 1, 0)
    }

  }

  return(CaseORControl)

}
