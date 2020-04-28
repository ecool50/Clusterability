run_monocle3Mod <- function(obj, clust_tech = "louvain",trueLabels, k.final){
  
  
  library(monocle)
  library(reticulate)
  import("louvain")
  
  # log
  data <- obj@Seurat$Scaled.Data
  col_info=data.frame(colnames(data))
  rownames(col_info)<- colnames(data)
  pd <- new("AnnotatedDataFrame", data = col_info)
  gene_info=data.frame(rownames(data))
  names(gene_info)="gene_short_name"
  rownames(gene_info)<- rownames(data)
  fd <- new("AnnotatedDataFrame", data = gene_info)
  cds <- monocle::newCellDataSet(cellData = data, phenoData = pd, featureData = fd)
  
  cds <- estimateSizeFactors(cds)
  # cds <- estimateDispersions(cds)
  
  cds <- monocle::preprocessCDS(cds, method = "PCA",  norm_method = "none", num_dim = k.final, relative_expr = FALSE, 
                                scaling = FALSE, pseudo_expr = 0) #Log normalize dataset
  
  cds <- monocle::reduceDimension(cds, reduction_method = "tSNE")
  
  cds <- monocle::clusterCells(cds, num_clusters = NULL, method = clust_tech, cores = 6) 
  
  log.clusters <- as.numeric(levels(cds$Cluster))[cds$Cluster]
  log.rand <- adjustedRandIndex(trueLabels, log.clusters)
  
  # multinom
  data <- obj@Scrna$Deviances
  col_info=data.frame(colnames(data))
  rownames(col_info)<- colnames(data)
  pd <- new("AnnotatedDataFrame", data = col_info)
  gene_info=data.frame(rownames(data))
  names(gene_info)="gene_short_name"
  rownames(gene_info)<- rownames(data)
  fd <- new("AnnotatedDataFrame", data = gene_info)
  cds <- monocle::newCellDataSet(cellData = data, phenoData = pd, featureData = fd)
  
  cds <- estimateSizeFactors(cds)
  # cds <- estimateDispersions(cds)
  
  cds <- monocle::preprocessCDS(cds, method = "PCA",  norm_method = "none", num_dim = k.final, relative_expr = FALSE, 
                                scaling = FALSE, pseudo_expr = 0) #Log normalize dataset
  
  cds <- monocle::reduceDimension(cds, reduction_method = "tSNE")
  
  cds <- monocle::clusterCells(cds, num_clusters = NULL, method = clust_tech, cores = 6) 
  
  multinom.clusters <- as.numeric(levels(cds$Cluster))[cds$Cluster]
  multinom.rand <- adjustedRandIndex(trueLabels, multinom.clusters)
  
  # negbinom
  data <- obj@SCTransform$Scaled.Data
  col_info=data.frame(colnames(data))
  rownames(col_info)<- colnames(data)
  pd <- new("AnnotatedDataFrame", data = col_info)
  gene_info=data.frame(rownames(data))
  names(gene_info)="gene_short_name"
  rownames(gene_info)<- rownames(data)
  fd <- new("AnnotatedDataFrame", data = gene_info)
  cds <- monocle::newCellDataSet(cellData = data, phenoData = pd, featureData = fd)
  
  cds <- estimateSizeFactors(cds)
  # cds <- estimateDispersions(cds)
  
  cds <- monocle::preprocessCDS(cds, method = "PCA",  norm_method = "none", num_dim = k.final, relative_expr = FALSE, 
                                scaling = FALSE, pseudo_expr = 0) #Log normalize dataset
  
  cds <- monocle::reduceDimension(cds, reduction_method = "tSNE")
  
  cds <- monocle::clusterCells(cds, num_clusters = NULL, method = clust_tech, cores = 6) 
  
  negbinom.clusters <- as.numeric(levels(cds$Cluster))[cds$Cluster]
  negbinom.rand <- adjustedRandIndex(trueLabels, negbinom.clusters)
  
  
  ensemble.mat <- cbind(log.clusters, multinom.clusters, negbinom.clusters)
  ensemble.labs <- diceR::LCA(ensemble.mat)
  ensemble.rand <- adjustedRandIndex(trueLabels, ensemble.labs)
  detach(package:monocle)
  return(list(Seurat = log.rand, Scrna = multinom.rand, SCT = negbinom.rand, Ensemble = ensemble.rand))
}