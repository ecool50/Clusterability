run_monocle3 <- function(data, dim_red = "tSNE", clust_tech = "louvain",trueLabels){
  
  
  library(monocle)
  library(reticulate)
  import("louvain")
  
  
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
  

  
  cds <- monocle::preprocessCDS(cds, method = "PCA",  norm_method = "log") #Log normalize dataset
  
  # if(use_dims_input=="TRUE"){
  #   cds <- monocle::reduceDimension(cds, reduction_method = dim_red, max_components=dims_input)
  # }
  
  
  cds <- monocle::reduceDimension(cds, reduction_method = "tSNE")

  
  
  cds <- monocle::clusterCells(cds, num_clusters = NULL, method = clust_tech, cores = 6) 
  
  clusters <- as.numeric(levels(cds$Cluster))[cds$Cluster]
  detach(package:monocle)
  return(list(ARI = adjustedRandIndex(trueLabels, clusters), clusters = clusters))
}