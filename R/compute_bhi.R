library(readr)
library(clusterProfiler)
library(GOSemSim)
library(org.Mm.eg.db)
source('bhi_functions.R')

# load the clusering table with each column gives the cluster ids of the samples in each row
partitions <- read.table("rep_average_breg.csv",sep=" ", header = TRUE, row.names=1)
genes <- readLines("gene1931.txt")
id_used <- intersect(which(genes!=""),which(duplicated(genes)==0))


######### BHI
mmGO <- godata('org.Mm.eg.db', ont="BP")

# get ENTREZID
genes_eid <- bitr(geneID = genes[id_used], fromType = 'SYMBOL',
                  toType = 'ENTREZID',OrgDb = 'org.Mm.eg.db',drop = FALSE)
genes_ENTREZID <- genes_eid$ENTREZID

# build sementic gene similarity matrix, very slow
g_matrix <- geneSimMatirx(genes_ENTREZID, mmGO, measure = "Resnik", 
                          drop = "IEA",combine = "max", verbose = FALSE)

n_g <- length(id_used)
scores <-  matrix(NA, nrow=n_g, ncol=n_g)
rownames(scores) <- genes_ENTREZID
colnames(scores) <- genes_ENTREZID
for (i in 1:n_g){
  for (j in 1:i){
    if(!is.na(genes_ENTREZID[i]) & !is.na(genes_ENTREZID[j]))
      scores[genes_ENTREZID[i],genes_ENTREZID[j]] <- g_matrix[genes_ENTREZID[i],genes_ENTREZID[j]]
  }
}
diag(scores) <- NA
g_matrix <- scores

# compute z_scores
n_time <- length(names(partitions))
zscores <- rep(0,n_time)
for(i in 1:n_time){
  clusters <- partitions[[i]][id_used]
  zscores[i] <- bhi_zscore(clusters,g_matrix,n_permute=1000)

}
View(t(zscores))

