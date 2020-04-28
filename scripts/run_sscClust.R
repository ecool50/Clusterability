run_sscClust<- function(data, dim_red = "pca", clust_tech = "hclust", trueLabels, ncomp = 25, k){
  
  library(sscClust)
  
  # if( preproc != "none" || (dim_red!="iCor" && dim_red!="PCA") || (use_dims_input!="TRUE" && use_dims_input!="internal") ||  (clust_tech!="k-means" && clust_tech!="hclust" && clust_tech!="SNN" && clust_tech!="ADPclust") || (clust!="set" && clust!="estimate")){
  #   stop("Provide correct parameters !")
  # }
  # 
  
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = data))
  rowData(sce)$feature_symbol <- rownames(data)
  
  #zinbwave dim red method doesnt work
  #SC3 clust method doesnt work 
  
  
  #clust=="set"
  # if( clust_tech=="kmeans" && clust=="set" && use_dims_input=="TRUE"){
  #   sscc <- sscClust::ssc.run(sce, assay.name="counts", method.reduction=dim_red, method.clust = clust_tech, k.batch=k_input, pca.npc=dims_input) 
  # }
  # 
  # else if( clust_tech=="kmeans" && clust=="set" && use_dims_input=="internal"){
  #   sscc <- sscClust::ssc.run(sce, assay.name="counts", method.reduction=dim_red, method.clust = clust_tech, k.batch=k_input) 
  # }
  # 
  # else if( clust_tech=="hclust" && clust=="set" && use_dims_input=="TRUE"){
  #   sscc <- sscClust::ssc.run(sce, assay.name="counts", method.reduction=dim_red, method.clust = clust_tech, k.batch=k_input, pca.npc=dims_input) 
  # }
  # 
  # else if( clust_tech=="hclust" && clust=="set" && use_dims_input=="internal"){
  #   sscc <- sscClust::ssc.run(sce, assay.name="counts", method.reduction=dim_red, method.clust = clust_tech, k.batch=k_input) 
  # }
  # 
  # else if( clust_tech=="adpclust" && clust=="set" && use_dims_input=="TRUE"){
  #   sscc <- sscClust::ssc.run(sce, assay.name="counts", method.reduction=dim_red, method.clust = clust_tech, k.batch=k_input, pca.npc=dims_input) 
  # }
  # 
  # else if( clust_tech=="adpclust" && clust=="set" && use_dims_input=="internal"){
  #   sscc <- sscClust::ssc.run(sce, assay.name="counts", method.reduction=dim_red, method.clust = clust_tech, k.batch=k_input) 
  # }
  
  
  
  #clust=="estimate"

  sscc <- sscClust::ssc.run(sce, assay.name="counts", method.reduction=dim_red, method.clust = clust_tech, ncore = 6, k.batch = k-5:k+5) 
  
  clusters <- as.numeric(factor(colData(sscc)[[1]])) #originally clusters are characters
  detach(package:sscClust)
  return(list(ARI = adjustedRandIndex(trueLabels, clusters), clusters = clusters))
} 