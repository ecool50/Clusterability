library(diptest)
library(dbscan)
library(MASS)
library(stats)
library(multiview)
library(mvc)
library(pcaReduce)
library(ascend)
library(TSCAN)
library(cidr)
library(flowCore)
library(FlowSOM)
library(monocle)
library(scran)
setwd("/home/elijah/Documents/skinny-dip/code/skinny-dip")

source("func.R")

data <- obj@Seurat$Sig.PCs.Data
pc <- data[, 1]
mode.intervals <- extractModalIntervals(sort(t2), 0.05, F)
resultLabels <- skinnyDipClusteringFullSpace(obj@Seurat$Umap.Data)


numNonNoiseClusters <- max(resultLabels)

# predictedLabels <- resultLabels

# run kmeans
# if(numNonNoiseClusters>1){
#   ## Run k-means if we have two clusters or more
#   clusterCenters <- matrix(nrow=0,ncol=ncol(data))
#   for(i in 1:numNonNoiseClusters){
#     dataObjectsInCluster <- data[predictedLabels==i,,drop=FALSE]
#     clusterCenters <- rbind(clusterCenters,colMeans(dataObjectsInCluster))
#   }
#   tryCatch({            
#     kmeansResult <- kmeans(data,clusterCenters)
#     predictedLabels <- kmeansResult$cluster        
#   },error=function(e){
#     if(grepl("empty cluster",as.character(e))){
#       print("We weren't able to initialize k-means with the centers, because there were one or more empty clusters that would have resulted. Reverting to the labels found by skinny-dip")
#     } else if(grepl("not distinct",as.character(e))){
#       print("We weren't able to initialize k-means with the centers, because the centers weren't distinct")
#     } else{
#       stop(e)
#     }
#   })        
# }

# plot(data[,1], data[,2], col=predictedLabels)


## Multiview
x1 <- obj@Seurat$Sig.PCs.Data
x2 <- obj@Scrna$Sig.PCs.Data
x3 <- t(obj@Seurat$Scaled.Data)
x4 <- t(obj@Scrna$Deviances)
x5 <- obj@Seurat$MDH$fitted
x6 <- obj@Scrna$MDH$fitted
x7 <- obj@Seurat$Tsne.Data
x8 <- obj@Scrna$Tsne.Data
x9 <- obj@Seurat$Umap.Data
x10 <- obj@Scrna$Umap.Data

t1 <- WGCNA::cor1(t(x1))
t2 <- WGCNA::cor1(t(x2))
t3 <- WGCNA::cor1(t(x3))
t4 <- WGCNA::cor1(t(x4))

# CIDR Preprocessing
sData <- scDataConstructor(obj@Data$cleaned, tagType = "raw")
sData <- determineDropoutCandidates(sData)
sData <- wThreshold(sData)
sData <- scDissim(sData, threads = 6)
sData <- scPCA(sData, plotPC = FALSE)
sData <- nPC(sData)

# Flowsom Preprocessing
dat <- logcounts(sce)
pca <- prcomp(t(dat), center = TRUE, scale. = FALSE)
pca <- pca$x[, seq_len(50), drop = FALSE]
ff <- flowFrame(exprs = pca)
fSOM <- FlowSOM::ReadInput(ff, compensate = FALSE, transform = FALSE, 
                           scale = FALSE, silent = TRUE)
fSOM <- FlowSOM::BuildSOM(fSOM, silent = TRUE, xdim = 15, 
                          ydim = 15)

# TSCAN 
dat <- obj@Data$counts[rowVars(obj@Data$counts) > 0, ]

# Monocle
cds <- convertTo(sce, type = "monocle")
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- reduceDimension(cds, max_components = 3, 
                       num_dim = 50,
                       reduction_method = "tSNE", verbose = TRUE)

# pcaReduce



numNonNoiseClusters <- 11
# res.dim <- multiview::mvmds(list(x3, x4), 2)
res.multi <- mvsc(list(x1, x2, x5, x6, x7, x8, x9, x10), 11)

res.cluster <- kmeansruns(res.multi$evectors, 
                          krange = numNonNoiseClusters, runs = 100, iter.max = 10 + round(log(nrow(x1))))$cluster
res.kmeans <- kmeansruns(data,  krange = numNonNoiseClusters, runs = 100,  iter.max = 10 + round(log(nrow(x1))))$cluster
res.pam <- pam(parallelDist(res.multi$evectors), diss = T, numNonNoiseClusters, metric = "euclidean")$clustering
res.hclust <- cutree(fastcluster::hclust(parallelDist(res.multi$evectors), method = "ward.D2"), numNonNoiseClusters)
cor.vals <- as.dist(1 - abs(WGCNA::cor1(t(res.multi$evectors))))
res.hclust.cor <- cutree(fastcluster::hclust(cor.vals, method = "complete"), k = numNonNoiseClusters)

res.kmeans <- kmeans(data, numNonNoiseClusters, nstart = 25)$cluster
res.hclust <- cutree(hclust(parallelDist(res.cluster$evectors), method = "ward.D2"), numNonNoiseClusters)
res.sc3 <- as.numeric(sce$sc3_9_clusters)
res.cidr <- scCluster(object = sData, nCluster = numNonNoiseClusters, nPC = sData@nPC, cMethod = "ward.D2")@clusters
res.flosom <- metaClustering_consensus(fSOM$map$codes, k = numNonNoiseClusters)
res.flosom <- res.flosom[fSOM$map$mapping[, 1]]
res.rtsneKmeans <- kmeans(obj@Seurat$Tsne.Data, centers = numNonNoiseClusters, nstart = 25)$cluster
res.umapKmeans <- kmeans(obj@Seurat$Umap.Data, centers = numNonNoiseClusters, nstart = 25)$cluster
res.tscan <- exprmclust(dat, clusternum = numNonNoiseClusters, modelNames = "VVV", reduce = TRUE)$clusterid
res.monocle <- clusterCells(cds, num_clusters = numNonNoiseClusters + 1, method = "densityPeak")$Cluster
