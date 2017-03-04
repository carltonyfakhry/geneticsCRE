# This functions joins the paths in Rels1 starting in Rels1_srcs to paths in Rels2
# starting in Rels2_srcs. The paths are joined by the target nodes in
# in Rels1, i.e Rels1_trgs, and the source nodes in Rels2, i.e Rels2_srcs. If
# path_length <= 3, then cummulation vectors are stored and returned to the user, otherwise
# they will be too big to fit in memory if path_length > 3.
ProcessPaths <- function(Rels1, Rels2,
                         Rels1_srcs, Rels1_trgs, joining_gene_sign, Rels2_srcs,
                         path_length, ValueTable, nCases, nControls, K, iterations,
                         CaseORControl, method, nthreads,
                         paths0, paths1, paths_res,
                         pos_path1, neg_path1, conflict_path1,
                         pos_path2, neg_path2, conflict_path2,
                         dest_path_pos1, dest_path_neg1, dest_path_conflict1){

  # Handle source nodes, which are source nodes in Rels1
  unique_Rels1_srcs <- sort(unique(Rels1_srcs))

  # Handle target nodes, which are source nodes in Rels2
  unique_Rels2_srcs <- sort(unique(Rels2_srcs))
  counts_Rels2_srcs <- rle(Rels2_srcs)$lengths
  cumsum_Rels2_srcs <- cumsum(counts_Rels2_srcs)
  locations_Rels2_srcs <- cumsum_Rels2_srcs - counts_Rels2_srcs # locations of the first occurrence of the ordered source nodes in Rels2

  # Get a mapping between uid and where it occurs in Rels2 and how many times it occurs contiguously
  # (p.s it could only occur contiguously since the source column in Rels2 is assumed to be ordered)
  uids_CountLoc <- geneticsCRE:::getMatchingList(unique_Rels2_srcs, counts_Rels2_srcs, locations_Rels2_srcs)

  inds <- which(!(unique(Rels1_trgs) %in% Rels2_srcs)) # account for the target nodes in Rels1 which are not present in the source nodes of Rels2
  for(i in unique(Rels1_trgs)[inds]){
    uids_CountLoc[[toString(i)]] = c(0,-1)
  }

  print(paste(length(unique_Rels1_srcs) , "source nodes for paths of length", path_length, "and their permutations will be processed!", sep = " "))

  lst <- geneticsCRE:::JoinIndices(Rels1_srcs, Rels1_trgs, uids_CountLoc, joining_gene_sign,
                                   ValueTable, nCases, nControls, K,
                                   iterations, CaseORControl, method, path_length, nthreads,
                                   paths0, paths1, paths_res,
                                   pos_path1, neg_path1, conflict_path1,
                                   pos_path2, neg_path2, conflict_path2,
                                   dest_path_pos1, dest_path_neg1, dest_path_conflict1)

  # lst <- geneticsCRE:::JoinIndicesMethod1(Rels1_srcs, Rels1_trgs, uids_CountLoc, joining_gene_sign,
  #                                         ValueTable, nCases, nControls, K,
  #                                         iterations, CaseORControl, method, path_length, nthreads,
  #                                         pos_path1, pos_path2, dest_path_pos1)

  # Get the paths and signedpaths
  inds <- lst$ids
  path <- c()
  signpath <- c()
  inds1 <- inds[,1]
  inds2 <- inds[,2]

  if(path_length == 1){

    Rels2_srcs2 <- c(Rels2_srcs, Rels2_srcs)
    path <- paste(Rels2_srcs2[inds2], sep = "")
    signpath <- c(paste(Rels2_srcs, " (+)", sep = ""), paste(Rels2_srcs, " (-)", sep = ""))
    signpath <- signpath[inds2]

  }else if(path_length == 2){

    inds3 <- inds2
    inds3[inds3 > nrow(Rels2)] <- inds3[inds3 > nrow(Rels2)] - nrow(Rels2)
    srcs <- Rels2$srcuid[inds3]
    trgs <- Rels2$trguid[inds3]
    path <- paste(srcs, trgs, sep = " -> ")

    dirs1 <- ifelse(inds2 > nrow(Rels2), "(-)", "(+)")

    dirs2 <- rep("(+)", length(inds2))
    sign1 <- Rels2$sign[inds3]
    dirs2[sign1 == 1 & dirs1 == "(-)"] <- "(-)"
    dirs2[sign1 == -1 & dirs1 == "(+)"] <- "(-)"

    signpath <- paste(paste(srcs, dirs1 , sep = " "),  paste(trgs, dirs2, sep = " "), sep = " -> ")

  }else if(path_length == 3){

    srcs <- Rels1$srcuid[inds1]
    trgs <- Rels1$trguid[inds1]
    inds3 <- inds2
    inds3[inds3 > nrow(Rels2)] <- inds3[inds3 > nrow(Rels2)] - nrow(Rels2)
    trgs2 <- Rels2$trguid[inds3]
    path <- paste(srcs, trgs, trgs2, sep = " -> ")

    dirs1 <- ifelse(inds2 > nrow(Rels1), "(-)", "(+)")

    dirs2 <- rep("(+)", length(inds1))
    sign1 <- Rels1$sign[inds1]
    dirs2[sign1 == 1 & dirs1 == "(-)"] <- "(-)"
    dirs2[sign1 == -1 & dirs1 == "(+)"] <- "(-)"

    dirs3 <- rep("(+)", length(inds1))
    sign2 <- Rels2$sign[inds3]
    dirs3[sign2 == 1 & dirs2 == "(-)"] <- "(-)"
    dirs3[sign2 == -1 & dirs2 == "(+)"] <- "(-)"

    signpath <- paste(paste(srcs, dirs1 , sep = " "),  paste(trgs, dirs2, sep = " "), paste(trgs2, dirs3, sep = " "), sep = " -> ")

  }else if(path_length == 4){

    srcs <- Rels1$srcuid[inds1]
    trgs <- Rels1$trguid[inds1]
    trgs2 <- Rels1$trguid2[inds1]
    inds3 <- inds2
    inds3[inds3 > nrow(Rels2)] <- inds3[inds3 > nrow(Rels2)] - nrow(Rels2)
    trgs3 <- Rels2$trguid[inds3]
    path <- paste(srcs, trgs, trgs2, trgs3, sep = " -> ")

    dirs1 <- ifelse(inds2 > nrow(Rels1), "(-)", "(+)")

    dirs2 <- rep("(+)", length(inds1))
    sign1 <- Rels1$sign[inds1]
    dirs2[sign1 == 1 & dirs1 == "(-)"] <- "(-)"
    dirs2[sign1 == -1 & dirs1 == "(+)"] <- "(-)"

    dirs3 <- rep("(+)", length(inds1))
    sign2 <- Rels1$sign2[inds1]
    dirs3[sign2 == 1 & dirs2 == "(-)"] <- "(-)"
    dirs3[sign2 == -1 & dirs2 == "(+)"] <- "(-)"

    dirs4 <- rep("(+)", length(inds1))
    sign3 <- Rels2$sign[inds3]
    dirs4[sign3 == 1 & dirs3 == "(-)"] <- "(-)"
    dirs4[sign3 == -1 & dirs3 == "(+)"] <- "(-)"

    signpath <- paste(paste(srcs, dirs1 , sep = " "),  paste(trgs, dirs2, sep = " "), paste(trgs2, dirs3, sep = " "), paste(trgs3, dirs4, sep = " "), sep = " -> ")

  }else if(path_length == 5){

    srcs <- Rels1$srcuid[inds1]
    trgs <- Rels1$trguid[inds1]
    trgs2 <- Rels1$trguid2[inds1]
    inds3 <- inds2
    inds3[inds3 > nrow(Rels2)] <- inds3[inds3 > nrow(Rels2)] - nrow(Rels2)
    trgs3 <- Rels2$trguid[inds3]
    trgs4 <- Rels2$trguid2[inds3]
    path <- paste(srcs, trgs, trgs2, trgs3, trgs4, sep = " -> ")

    dirs1 <- ifelse(inds2 > nrow(Rels1), "(-)", "(+)")

    dirs2 <- rep("(+)", length(inds1))
    sign1 <- Rels1$sign[inds1]
    dirs2[sign1 == 1 & dirs1 == "(-)"] <- "(-)"
    dirs2[sign1 == -1 & dirs1 == "(+)"] <- "(-)"

    dirs3 <- rep("(+)", length(inds1))
    sign2 <- Rels1$sign2[inds1]
    dirs3[sign2 == 1 & dirs2 == "(-)"] <- "(-)"
    dirs3[sign2 == -1 & dirs2 == "(+)"] <- "(-)"

    dirs4 <- rep("(+)", length(inds3))
    sign3 <- Rels2$sign[inds3]
    dirs4[sign3 == 1 & dirs3 == "(-)"] <- "(-)"
    dirs4[sign3 == -1 & dirs3 == "(+)"] <- "(-)"

    dirs5 <- rep("(+)", length(inds3))
    sign4 <- Rels2$sign2[inds3]
    dirs5[sign4 == 1 & dirs4 == "(-)"] <- "(-)"
    dirs5[sign4 == -1 & dirs4 == "(+)"] <- "(-)"

    signpath <- paste(paste(srcs, dirs1 , sep = " "),  paste(trgs, dirs2, sep = " "), paste(trgs2, dirs3, sep = " "), paste(trgs3, dirs4, sep = " "), paste(trgs4, dirs5, sep = " "), sep = " -> ")

  }


  lst2 <- list()
  lst2[["scores"]] <- lst$scores
  lst2[["TestScores"]] <- lst$TestScores
  lst2[["path"]] <- path
  lst2[["signpath"]] <- signpath
  if(path_length <= 3){
    lst2[["cummulation_vectors"]] <- lst[["cummulation_vectors"]]
  }
  if(path_length == 3 | path_length == 4){
    inds2 <- inds[,2]
    inds2[inds2 > nrow(Rels2)] <- inds2[inds2 > nrow(Rels2)] - nrow(Rels2)
    inds[,2] <- inds2
    lst2[["inds"]] <- inds
  }

  return(lst2)

}
