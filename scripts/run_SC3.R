run_SC3 <- function(data, ncells = 500, trueLabels){
  library(SC3)
  sce <- SingleCellExperiment(
    assays = list(counts = data
    )
  )
  sce <- scater::normalize(sce)
  sce <- SC3::sc3_estimate_k(sce)
  sc3_num_clust <- sce@metadata$sc3$k_estimation
  if(sc3_num_clust == 1){
    return(list(ARI = adjustedRandIndex(trueLabels, rep(1,ncol(data))), clusters = rep(1,ncol(data))))
  } 
  
  # define feature names in feature_symbol column
  rowData(sce)$feature_symbol <- rownames(sce)
  sce <- SC3::sc3(sce, ks =sc3_num_clust,  biology = F, svm_num_cells = ncells, gene_filter = T, n_cores = 6)
  sce <- sc3_run_svm(sce, ks = sc3_num_clust)
  res <- as.numeric(colData(sce)[, paste0("sc3_", sc3_num_clust, "_clusters")])
  # res[is.na(res)] <- -1
  return(list(ARI = adjustedRandIndex(trueLabels, res), clusters = res))
}