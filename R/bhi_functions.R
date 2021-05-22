##### functions

# modified from mgeneSim in GOSemSim
geneSimMatirx <- function (genes, semData, measure="Resnik", drop="IEA", combine="max", verbose=FALSE) {
  genes <- unique(as.character(genes))
  n <- length(genes)
  scores <- matrix(NA, nrow=n, ncol=n)
  rownames(scores) <- genes
  colnames(scores) <- genes
  
  gos <- lapply(genes, gene2GO, godata=semData, dropCodes=drop)
  uniqueGO <-  unique(unlist(gos))
  go_matrix <- mgoSim(uniqueGO, uniqueGO, semData, measure = measure, combine = NULL)
  if (verbose) {
    cnt <- 1
    pb <- txtProgressBar(min=0, max=sum(1:n), style=3)
  }
  for (i in seq_along(genes)) {
    for (j in seq_len(i)){
      if (verbose) {
        setTxtProgressBar(pb, cnt)
        cnt <- cnt + 1
      }
      scores[i,j] <- combineScores(go_matrix[gos[[i]], gos[[j]]], combine = combine)
      scores[j,i] <- scores[i,j]
    }
  }
  if (verbose)
    close(pb)
  #  removeRowNA <- apply(!is.na(scores), 1, sum)>0
  #  removeColNA <- apply(!is.na(scores), 2, sum)>0
  return(scores)
}

# copied from GOSemSim
gene2GO <- function(gene, godata, dropCodes) {
  goAnno <- godata@geneAnno
  if (! "EVIDENCE" %in% colnames(goAnno)) {
    warning("Evidence codes not found, 'drop' parameter will be ignored...")
  } else {
    goAnno <- goAnno[!goAnno$EVIDENCE %in% dropCodes,]
  }
  go <- as.character(unique(goAnno[goAnno[,1] == gene, "GO"]))
  go[!is.na(go)]
}

# Compute the sementic BHI for one partition
# @param clusters: clustering ids
# @param bhi_matrix: may contain NA as there are unmapped genes or genes with no GO terms.
cluster_BHI <- function (clusters, bhi_matrix){
  n_c <- length(unique(clusters))
  cluster_BHI <- matrix(NA,nrow = n_c)
  for(i in 1:n_c){
    idx <- which(clusters==i)
    cluster_BHI[i] <- mean(bhi_matrix[idx,idx],na.rm=TRUE)
  }
  return(mean(cluster_BHI,na.rm=TRUE))
}

# Compute the BHI z_score for one partition
# @param clusters: clustering ids
# @param bhi_matrix: may contain NA as there are unmapped genes or genes with no GO terms.
# @param n_permute: number of random clusterings
bhi_zscore <-function(clusters,bhi_matrix,n_permute = 500){
  random_bhi <- rep(0,n_permute)
  cluster_bhi <- cluster_BHI(clusters,bhi_matrix)
  for(i in 1:n_permute){
    cids = sample(clusters)
    random_bhi[i] <- cluster_BHI(cids,bhi_matrix)
  }
  return( (cluster_bhi - mean(random_bhi))/sd(random_bhi) )
}