#' This function runs the causal engine for a given phenotype dataset.
#'
#' @description This function runs the causal engine for a given phenotype dataset.
#'
#' @usage GetBestPaths(dataset, nCases, nControls, method = 1,
#'        threshold_percent = 0.05, K = 10, pathLength = 5, iterations = 100,
#'        strataF = NA)
#'
#' @param dataset Filename of the phenotype data. First column must correspond
#'                to the gene symbols; all other columns correspond to the patients.
#'                First column must be named symbols and must contrain strings corresponding
#'                to the symbols of the genes. All other columns either contain a 0 if the
#'                gene is not present in the patient, and 1 or 2 if the gene is present
#'                in the patient.
#' @param nCases The number of cases (Must be greater than 1).
#' @param nControls The number of controls (Must be greater than 1).
#' @param path The path where the package will write and read from files
#'              created by the package.
#' @param method The method to be invoked. The value 1 calls the unsigned method
#'               while the value 2 calls the signed method.
#' @param threshold_percent All genes which occur in less than emph{threshold_percent}
#'                          of patients are kept for processing. The default
#'                          value is 0.05.
#' @param K The top K paths to be returned of each path length. The default
#'          value is 10.
#' @param pathLength The maximum path length for which the top K paths
#'                   will be computed. The default value is 5 (The maximum
#'                   value for \emph{pathLength} is 5).
#' @param iterations The number of iterations to be used in computing the
#'                   p-values. The default value is 100 (The minimum acceptable
#'                   value is 1).
#' @param strataF Filename for the strata file. \emph{strataF} must be a file
#'                of two columns, first column named \emph{subjid }contains
#'                the patient names similar to the dataset columns, and the
#'                second column named \emph{stratum} must contain the integer
#'                corresponding to the stratum.
#' @param nthreads The number of threads to be used in the parallel region
#'                 of the code. Default value is \emph{nthreads = NA}, in
#'                 which case \emph{nthreads} will be set to the maximum number of
#'                 threads that can be created. If \emph{nthreads} is greater
#'                 than the maximum number of threads that can be created, then
#'                 \emph{nthreads} will default back to the maximum number of threads
#'                 that can be created.
#'
#' @return This function returns a list with the following items:
#'         \item{SignedPaths}{The top K signed paths for each length.}
#'         \item{UserKpaths}{The top K paths for each length.}
#'         \item{UserLengths}{The lengths of each path.}
#'         \item{UserScores}{The scores of each path.}
#'         \item{UserPvalues}{The p-values of each path.}
#'
#' @author Carl Tony Fakhry
#'
#' @references Franceschini, A (2013). STRING v9.1: protein-protein interaction networks, with increased coverage
#'             and integration. In:'Nucleic Acids Res. 2013 Jan;41(Database issue):D808-15. doi: 10.1093/nar/gks1094.
#'             Epub 2012 Nov 29'.
#'
#'@export
GetBestPaths <- function(dataset, nCases, nControls, path = ".", method = 1, threshold_percent = 0.05, K = 10, pathLength = 5, iterations = 100, strataF = NA, nthreads = NA){

  # Check the input parameters
  geneticsCRE:::check_input(nCases, nControls, method, threshold_percent, K, pathLength, iterations, nthreads)
  if(is.na(nthreads)) nthreads <- -1
  path <- normalizePath(path)
  if(substr(path, nchar(path),1) != "/") path <- paste(path, "/", sep = "")

  # Precompute values table
  print("Computing Values Table...")
  ValueTable <- geneticsCRE:::getValuesTable(nCases, nControls)

  # Preprocessing Phenotype dataset, all genes which occur in more
  # than threshold_percent of patients will be removed from the
  # dataset, also genes which do not occur in any patients will be removed
  print("Preprocessing Phenotype dataset...")
  dataset <- geneticsCRE:::PreprocessTable(dataset, threshold_percent, nCases, nControls)
  genes <- dataset$genes
  data <- dataset$data
  genes_data <- data.frame(genes = genes, data)
  patients <- dataset$patients
  rm(dataset)

  # Read in Stringdb (Already parsed)
  print("Processing network and dataset...")
  KB <- getStringKB()
  Rels <- KB$Rels
  Ents <- KB$Ents
  rm(KB)

  # # max uid is 8303, so GENE1 has uid 8304 UP TO 8308
  # new_uids <- (max(Ents$uid)+1):(max(Ents$uid)+5)
  # new_ents <- data.frame(uid = new_uids, id = 1:5, symbol = c("GENE1","GENE2","GENE3","GENE4","GENE5"),ensembleid=rep("",5))
  # Ents <- rbind(Ents, new_ents)
  #
  # #Experiment 1: Show method 2 and 3 work
  # Rels <- rbind(Rels, c(8304,8305,-1))
  # Rels <- rbind(Rels, c(8305,8306,-1))
  # Rels <- rbind(Rels, c(8306,8307,-1))
  # Rels <- rbind(Rels, c(8307,8308,-1))
  #
  # # The significant path should be GENE1 (+) -> GENE2 (-) -> GENE3 (+) -> GENE4 (-) -> GENE5 (+)  for method 3 and its unsigned version for method 2
  #
  # df = matrix(0, nrow = 5, ncol = ncol(data))
  #
  # # Add 5 ones for each protein
  # df[1,2:6] = 1 ; df[2,(1537+1):(1537+5)] = 1 ; df[3,401:405] = 1 ; df[4,(1537+201):(1537+205)] = 1 ; df[5,801:805] = 1
  #
  # # Add 5 more ones for each protein
  # df[1,7:11] = 1 ; df[2,(1537+6):(1537+10)] = 1 ; df[3,406:410] = 1 ; df[4,(1537+206):(1537+210)] = 1 ; df[5,806:810] = 1
  #
  # # Add 10 more ones for each protein
  # df[1,12:21] = 1 ; df[2,(1537+11):(1537+20)] = 1 ; df[3,411:420] = 1 ; df[4,(1537+211):(1537+220)] = 1 ; df[5,811:820] = 1
  #
  # # Add 20 more ones for each protein
  # df[1,22:41] = 1 ; df[2,(1537+21):(1537+40)] = 1 ; df[3,421:440] = 1 ; df[4,(1537+221):(1537+240)] = 1 ; df[5,821:840] = 1
  #
  # # Add 20 more ones for each protein
  # df[1,42:61] = 1 ; df[2,(1537+41):(1537+60)] = 1 ; df[3,441:460] = 1 ; df[4,(1537+241):(1537+260)] = 1 ; df[5,841:860] = 1
  #
  # # Add 40 more ones for each protein
  # df[1,62:101] = 1 ; df[2,(1537+61):(1537+100)] = 1 ; df[3,461:500] = 1 ; df[4,(1537+261):(1537+300)] = 1 ; df[5,861:900] = 1
  #
  # # Add 40 more ones for each protein
  # df[1,102:141] = 1 ; df[2,(1537+101):(1537+140)] = 1 ; df[3,501:540] = 1 ; df[4,(1537+301):(1537+340)] = 1 ; df[5,901:940] = 1


  # # Experiment 2 flipped version of Experiment1
  # # Add 5 ones for each protein
  # df[1,(1537+2):(1537+6)] = 1 ; df[2,1:5] = 1 ; df[3,(1537+401):(1537+405)] = 1 ; df[4,201:205] = 1 ; df[5,(1537+801):(1537+805)] = 1
  #
  # # Add 5 more ones for each protein
  # df[1,(1537+7):(1537+11)] = 1 ; df[2,6:10] = 1 ; df[3,(1537+406):(1537+410)] = 1 ; df[4,206:210] = 1 ; df[5,(1537+806):(1537+810)] = 1
  #
  # # Add 10 more ones for each protein
  # df[1,(1537+12):(1537+21)] = 1 ; df[2,11:20] = 1 ; df[3,(1537+411):(1537+420)] = 1 ; df[4,211:220] = 1 ; df[5,(1537+811):(1537+820)] = 1
  # #
  # # # Add 20 more ones for each protein
  # df[1,(1537+22):(1537+41)] = 1 ; df[2,21:40] = 1 ; df[3,(1537+421):(1537+440)] = 1 ; df[4,221:240] = 1 ; df[5,(1537+821):(1537+840)] = 1
  #
  # # Add 20 more ones for each protein
  # df[1,(1537+42):(1537+61)] = 1 ; df[2,41:60] = 1 ; df[3,(1537+441):(1537+460)] = 1 ; df[4,241:260] = 1 ; df[5,(1537+841):(1537+860)] = 1
  #
  # # # Add 40 more ones for each protein
  # df[1,(1537+62):(1537+101)] = 1 ; df[2,61:100] = 1 ; df[3,(1537+461):(1537+500)] = 1 ; df[4,261:300] = 1 ; df[5,(1537+861):(1537+900)] = 1
  #
  # # Add 40 more ones for each protein
  # df[1,(1537+102):(1537+141)] = 1 ; df[2,101:140] = 1 ; df[3,(1537+501):(1537+540)] = 1 ; df[4,301:340] = 1 ; df[5,(1537+901):(1537+940)] = 1


  # Experiment 2 Show method 3 can work when method 2 can't
  # Add 5 ones for each protein
  # df[2,(1537+1):(1537+5)] = 1 ; df[4,(1537+201):(1537+205)] = 1
  #
  # # # Add 5 more ones for each protein
  # df[2,(1537+6):(1537+10)] = 1 ; df[4,(1537+206):(1537+210)] = 1
  #
  # # # Add 10 more ones for each protein
  # df[2,(1537+11):(1537+20)] = 1 ; df[4,(1537+211):(1537+220)] = 1
  #
  # # Add 20 more ones for each protein
  # df[2,(1537+21):(1537+40)] = 1 ; df[4,(1537+221):(1537+240)] = 1
  #
  # # Add 20 more ones for each protein
  # df[2,(1537+41):(1537+60)] = 1 ; df[4,(1537+241):(1537+260)] = 1
  #
  # # Add 40 more ones for each protein
  # df[2,(1537+61):(1537+100)] = 1 ; df[4,(1537+261):(1537+300)] = 1
  #
  # # # Add 40 more ones for each protein
  # df[2,(1537+101):(1537+140)] = 1 ; df[4,(1537+301):(1537+340)] = 1

  # Statistical Significance for partial paths
  # df[1,2:141] = 1
  # df[2,(1537+1):(1537+140)] = 1
  # df[3,401:540] = 1
  # df[4,(1537+201):(1537+340)] = 1
  # df[5,801:940] = 1

  # data = rbind(data, df)
  # genes = c(genes, "GENE1","GENE2","GENE3","GENE4","GENE5")
  # genes_data <- data.frame(genes = genes, data)

  # Keep gene symbols in Ents which are present in the dataset
  Ents <- Ents[Ents$symbol != "-1",]
  Ents <- Ents[!duplicated(Ents$symbol),]
  Ents <- Ents[Ents$symbol %in% genes_data$genes,]
  Ents <- Ents[order(Ents$uid),]

  # Keep gene symbols in data which are present in Ents
  genes_data <- genes_data[genes_data$genes %in% Ents$symbol,]

  # Remove edges in Stringdb which have genes not present
  # in Ents
  Rels <- Rels[,c(1,2,3)]
  names(Rels) <- c("srcuid", "trguid", "sign")
  Rels <- unique(Rels)
  Rels <- Rels[Rels$srcuid %in% Ents$uid,]
  Rels2 <- Rels   # Rels2 could contain target genes which are not present in Ents (and dataset)
  Rels <- Rels[Rels$trguid %in% Ents$uid,] # Rels only contains both source and target genes present in Ents (and dataset)
  Rels <- Rels[Rels$srcuid != Rels$trguid,]

  # genes2, data2 and geneuids2 have the same order as Ents2$symbol
  # Keep only genes which are still present in source nodes of Network
  leftuids2 <- unique(Rels2$srcuid)
  Ents2 <- Ents[Ents$uid %in% leftuids2,]
  Ents2 <- Ents2[order(Ents2$uid),]
  genes_data2 <- genes_data[match(Ents2$symbol, genes_data$genes),]
  data2 <- as.matrix(genes_data2[,2:ncol(genes_data2),drop=FALSE])
  genes_data2 <- cbind(uid = Ents2$uid, genes_data2)

  # genes, data and geneuids have the same order as Ents$symbol
  # Keep only genes which are still present in Network
  leftuids <- unique(c(Rels$srcuid, Rels$trguid))
  Ents <- Ents[Ents$uid %in% leftuids,]
  Ents <- Ents[order(Ents$uid),]
  genes_data <- genes_data[match(Ents$symbol, genes_data$genes),]
  data <- as.matrix(genes_data[,2:ncol(genes_data),drop=FALSE])
  genes_data <- cbind(uid = Ents$uid, genes_data)

  stratagroups <- NA
  stratanumbers <- NA

  # Handling Strata File
  if(length(strataF) == 1 & !is.na(strataF)){
    strata <- read.table(strataF, header = T,stringsAsFactors = F)
    if(names(strata)[1] != "subjid") stop("First column name in the strata file should be named subjid!")
    if(names(strata)[2] != "stratum") stop("Second column name in the strata file should be named stratum!")
    if(!all(strata[,1] %in% patients)) stop("All patient names in the strata file should be present in the column names of dataset!")
    if(!all(patients %in% strata[,1])) stop("All patient names in the dataset should be present in the strata file!")
    inds <- match(patients, strata[,1])
    strata <- strata[inds,] # strata patients now have the same order as the column names in the dataset
    stratanumbers <- as.integer(as.vector(strata[,2]))
    stratagroups <- as.vector(unique(stratanumbers))
  }

  # Precompute the permuted columns to be used in producing randomized paths for computing
  # the p-values of the paths
  RandIndicesMat <- geneticsCRE:::getRandIndicesMat(ncol(data), iterations, strataF, nCases, nControls, stratagroups, stratanumbers)
  CaseORControl <- geneticsCRE:::getCaseORControl(RandIndicesMat, nCases, nControls)

  UserKpaths <- c()
  UserKsignpaths <- c()
  UserSymbPaths <- c()
  UserSymbSignPaths <- c()
  UserScores <- c()
  UserLengths <- c()
  UserPvalues <- c()

  # Create the relations data frame for the input data. This data frame will be used in the joining
  # of the paths of length 1 and 2.
  Rels_data <- data.frame(srcuid = Ents$uid)
  Rels_data2 <- data.frame(srcuid = Ents2$uid)

  # Create the relations data frame for paths of length 3
  Rels <- Rels[order(Rels$srcuid, Rels$trguid),]
  Rels2 <- Rels
  names(Rels2) <- c("trguid", "trguid2", "sign2")
  Rels3 <- dplyr::left_join(Rels, Rels2, by = "trguid")
  Rels3 <- Rels3 %>% filter(complete.cases(.))
  Rels3 <- Rels3[order(Rels3$srcuid, Rels3$trguid, Rels3$trguid2),]
  rm(Rels2)

  # Get the sign of the third gene in a path of length 3 assuming the first
  # gene in the path has a positive sign
  third_gene_sign <- rep(1, nrow(Rels3))
  third_gene_sign[Rels3$sign == -1 & Rels3$sign2 == 1] <- -1
  third_gene_sign[Rels3$sign == 1 & Rels3$sign2 == -1] <- -1

  dest_path_pos1 <- paste(path, "geneticsCRE_joined_pos1", sep = "")
  dest_path_neg1 <- paste(path, "geneticsCRE_joined_neg1", sep = "")
  dest_path_conflict1 <- paste(path, "geneticsCRE_joined_conflict1", sep = "")
  dest_path_pos2 <- paste(path, "geneticsCRE_joined_pos2", sep = "")
  dest_path_neg2 <- paste(path, "geneticsCRE_joined_neg2", sep = "")
  dest_path_conflict2 <- paste(path, "geneticsCRE_joined_conflict2", sep = "")
  dest_path_pos3 <- paste(path, "geneticsCRE_joined_pos3", sep = "")
  dest_path_neg3 <- paste(path, "geneticsCRE_joined_neg3", sep = "")
  dest_path_conflict3 <- paste(path, "geneticsCRE_joined_conflict3", sep = "")

  timeStart <- as.numeric(Sys.time()) * 1000

  # Processing the paths
  for(path_length in 1:pathLength){

    if(path_length == 1){

      print("Processing paths of length 1...")

      path_data_file <- paste(path, "data_geneticsCRE", sep = "")
      geneticsCRE:::parsePaths(data, nCases, nControls, path_data_file)
      rm(data)

      lst1 <- geneticsCRE:::ProcessPaths(Rels_data, Rels_data, Rels_data$srcuid, Rels_data$srcuid, rep(1,length(Ents$uid)), Rels_data$srcuid,
                                         1, ValueTable, nCases, nControls, K, iterations, CaseORControl, method, nthreads, "", "", "", path_data_file, "", "",dest_path_pos1,
                                         dest_path_neg1, dest_path_conflict1)
      file.remove(path_data_file)


      path_data2_file <- paste(path, "data2_geneticsCRE", sep = "")
      geneticsCRE:::parsePaths(data2, nCases, nControls, path_data2_file)
      rm(data2)

      lst <- geneticsCRE:::ProcessPaths(Rels_data2, Rels_data2, Rels_data2$srcuid, Rels_data2$srcuid, rep(1,length(Ents2$uid)), Rels_data2$srcuid,
                                        1, ValueTable, nCases, nControls, K, iterations, CaseORControl, method, nthreads, "", "", "", path_data2_file, "", "", "",
                                        "", "")
      file.remove(path_data2_file)

    }else if(path_length == 2){

      print("Processing paths of length 2...")

      path_data3_file <- paste(path, "data3_geneticsCRE", sep = "")
      data3 <- geneticsCRE:::matchData(genes_data, Rels$trguid)
      geneticsCRE:::parsePaths(data3, nCases, nControls, path_data3_file)
      rm(data3)

      lst2 <- geneticsCRE:::ProcessPaths(Rels_data, Rels, Rels_data$srcuid, Rels_data$srcuid, Rels$sign, Rels$srcuid,
                                          2, ValueTable, nCases, nControls, K, iterations, CaseORControl, method, nthreads, dest_path_pos1,
                                         dest_path_neg1, dest_path_conflict1, path_data3_file, "", "", dest_path_pos2,
                                         dest_path_neg2, dest_path_conflict2)
      lst <- lst2
      file.remove(path_data3_file)

    }else if(path_length == 3){

      print("Processing paths of length 3...")

      path_data4_file <- paste(path, "data4_geneticsCRE", sep = "")
      data4 <- geneticsCRE:::matchData(genes_data, Rels$trguid)
      geneticsCRE:::parsePaths(data4, nCases, nControls, path_data4_file)
      rm(data4)

      lst3 <- geneticsCRE:::ProcessPaths(Rels, Rels, Rels$srcuid, Rels$trguid, Rels$sign, Rels$srcuid,
                                         3, ValueTable, nCases, nControls, K, iterations, CaseORControl, method, nthreads,
                                         dest_path_pos2, dest_path_neg2, dest_path_conflict2, path_data4_file, "", "",
                                         dest_path_pos3, dest_path_neg3, dest_path_conflict3)

      lst <- lst3
      file.remove(path_data4_file)

    }
    else if(path_length == 4){

      print("Processing paths of length 4...")


      lst4 <- geneticsCRE:::ProcessPaths(Rels3, Rels, Rels3$srcuid, Rels3$trguid2, third_gene_sign, Rels$srcuid,
                                         4, ValueTable, nCases, nControls, K, iterations,
                                         CaseORControl, method, nthreads,
                                         dest_path_pos3, dest_path_neg3, dest_path_conflict3,
                                         dest_path_pos2, dest_path_neg2, dest_path_conflict2,
                                         "", "", "")
      lst <- lst4

    }
    else if(path_length == 5){


      print("Processing paths of length 5...")

      lst5 <- geneticsCRE:::ProcessPaths(Rels3, Rels3, Rels3$srcuid, Rels3$trguid2, third_gene_sign, Rels3$srcuid,
                                         5, ValueTable, nCases, nControls, K, iterations,
                                         CaseORControl, method, nthreads,
                                         dest_path_pos3, dest_path_neg3, dest_path_conflict3,
                                         dest_path_pos3, dest_path_neg3, dest_path_conflict3,
                                         "", "", "")

      lst <- lst5

    }

    UserScores <- c(UserScores, lst$scores)
    UserKpaths <- c(UserKpaths, lst$path)
    if(path_length == 1){
      UserSymbPaths <- c(UserSymbPaths, geneticsCRE:::Uid2Symbol(Ents2, lst$path))
    }else{
      UserSymbPaths <- c(UserSymbPaths, geneticsCRE:::Uid2Symbol(Ents, lst$path))
    }
    UserKsignpaths <- c(UserKsignpaths, lst$signpath)
    if(path_length == 1){
      UserSymbSignPaths <- c(UserSymbSignPaths, geneticsCRE:::SignUid2SignSymbol(Ents2, lst$signpath))
    }else{
      UserSymbSignPaths <- c(UserSymbSignPaths, geneticsCRE:::SignUid2SignSymbol(Ents, lst$signpath))
    }
    TestScores <- lst$TestScores

    Pvalues <- sapply(lst$scores, function(x, y) {length(which(TestScores >= x))/length(TestScores)}, y = TestScores)
    UserPvalues <- c(UserPvalues, Pvalues)
    UserLengths <- c(UserLengths, rep(path_length, length(Pvalues)))

  }

  print(sprintf("total time: %f", as.numeric(Sys.time()) * 1000 - timeStart, digits=15))

  suppressWarnings(file.remove(dest_path_pos1, dest_path_neg1, dest_path_conflict1))
  suppressWarnings(file.remove(dest_path_pos2, dest_path_neg2, dest_path_conflict2))
  suppressWarnings(file.remove(dest_path_pos3, dest_path_neg3, dest_path_conflict3))

  BestPaths <- data.frame(UserSymbSignPaths, UserSymbPaths, UserLengths, UserScores, UserPvalues, stringsAsFactors = F)
  BestPaths <- BestPaths[order(BestPaths$UserPvalues),]
  names(BestPaths) <- c("SignedPaths", "Paths", "Lengths", "Scores", "Pvalues")

  return(BestPaths)

}
