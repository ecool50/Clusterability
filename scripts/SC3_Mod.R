run_SC3Mod <- function(obj, ncells = 500, trueLabels){
  #seurat
  sce <- SingleCellExperiment(
    assays = list(counts = obj@Seurat$Scaled.Data, logcounts = obj@Seurat$Scaled.Data
    )
  )
  sce <- SC3::sc3_estimate_k(sce)
  sc3_num_clust <- sce@metadata$sc3$k_estimation
  
  # define feature names in feature_symbol column
  rowData(sce)$feature_symbol <- rownames(sce)
  sce <- SC3::sc3(sce, ks =sc3_num_clust ,svm_num_cells = ncells,  biology = F, gene_filter = F, n_cores = 6)
  sce <- SC3::sc3_run_svm(sce, ks = sc3_num_clust)
  res.seurat <- as.numeric(colData(sce)[, paste0("sc3_", sc3_num_clust, "_clusters")])
  seurat.rand <- adjustedRandIndex(trueLabels, res.seurat)
  gc()
  # Scrna
  sce <- SingleCellExperiment(
    assays = list(counts = obj@Scrna$Deviances, logcounts = obj@Scrna$Deviances
    )
  )
  sce <- SC3::sc3_estimate_k(sce)
  sc3_num_clust <- sce@metadata$sc3$k_estimation
  
  # define feature names in feature_symbol column
  rowData(sce)$feature_symbol <- rownames(sce)
  sce <- SC3::sc3(sce, ks =sc3_num_clust ,svm_num_cells = ncells,  biology = F, gene_filter = F, n_cores = 6)
  sce <- SC3::sc3_run_svm(sce, ks = sc3_num_clust)
  res.scrna <- as.numeric(colData(sce)[, paste0("sc3_", sc3_num_clust, "_clusters")])
  scrna.rand <- adjustedRandIndex(trueLabels, res.scrna)
  gc()
  #SCT
  sce <- SingleCellExperiment(
    assays = list(counts = obj@SCTransform$Scaled.Data, logcounts = obj@SCTransform$Scaled.Data
    )
  )
  sce <- SC3::sc3_estimate_k(sce)
  sc3_num_clust <- sce@metadata$sc3$k_estimation
  
  # define feature names in feature_symbol column
  rowData(sce)$feature_symbol <- rownames(sce)
  sce <- SC3::sc3(sce, ks =sc3_num_clust ,svm_num_cells = ncells,  biology = F, gene_filter = F, n_cores = 6)
  sce <- SC3::sc3_run_svm(sce, ks = sc3_num_clust)
  res.sct<- as.numeric(colData(sce)[, paste0("sc3_", sc3_num_clust, "_clusters")])

  sct.rand <- adjustedRandIndex(trueLabels, res.sct)
  
  ensemble.mat <- cbind(res.seurat, res.scrna, res.sct)
  ensemble <- apply(ensemble.mat, MARGIN = 2, as.cl_partition)
  cons <- cl_consensus(as.cl_ensemble(ensemble), method = "SE", control = list(nruns = 50))
  ensemble.labs <- c(cl_class_ids(cons))
  ensemble.rand <- adjustedRandIndex(trueLabels, ensemble.labs)
  gc()
  return(list(Seurat = seurat.rand, Scrna = scrna.rand, SCT = sct.rand, Ensemble = ensemble.rand))
}