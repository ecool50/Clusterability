Run_sscClustMod <- function(obj, ncomp.seurat, ncomp.scrna, ncomp.sct, trueLabels, dim_red = "pca", clust_tech = "hclust", k, truelabels){
  #seurat
  sce <- SingleCellExperiment(
    assays = list(counts = obj@Seurat$Scaled.Data, logcounts = obj@Seurat$Scaled.Data
    )
  )
  rowData(sce)$feature_symbol <- rownames(sce)
  sscc <- sscClust::ssc.run(sce, assay.name="counts", method.reduction=dim_red, method.clust = clust_tech, pca.npc=k.seurat, ncore = 6,
                            k.batch = k-2:k+2) 
  
  res.seurat <- as.numeric(factor(colData(sscc)[[1]])) #originally clusters are characters
  seurat.rand <- adjustedRandIndex(truelabels, res.seurat)
  
  # Scrna
  sce <- SingleCellExperiment(
    assays = list(counts = obj@Scrna$Deviances, logcounts = obj@Scrna$Deviances
    )
  )
  rowData(sce)$feature_symbol <- rownames(sce)
  sscc <- sscClust::ssc.run(sce, assay.name="counts", method.reduction=dim_red, method.clust = clust_tech, pca.npc=k.scrna, ncore = 6,
                            k.batch = k-2:k+2) 
  
  res.scrna <- as.numeric(factor(colData(sscc)[[1]])) #originally clusters are characters
  scrna.rand <- adjustedRandIndex(truelabels, res.scrna)
  
  #SCT
  sce <- SingleCellExperiment(
    assays = list(counts = obj@SCTransform$Scaled.Data, logcounts = obj@SCTransform$Scaled.Data
    )
  )
  rowData(sce)$feature_symbol <- rownames(sce)
  sscc <- sscClust::ssc.run(sce, assay.name="counts", method.reduction=dim_red, method.clust = clust_tech, pca.npc=k.sct, ncore = 6,
                            k.batch = k-2:k+2) 
  
  res.sct <- as.numeric(factor(colData(sscc)[[1]])) #originally clusters are characters
  sct.rand <- adjustedRandIndex(truelabels, res.sct)
  
  # compute ensemble
  ensemble.mat <- cbind(res.seurat, res.scrna, res.sct)
  ensemble.labs <- diceR::LCA(ensemble.mat)
  ensemble.rand <- adjustedRandIndex(truelabels, ensemble.labs)
  return(list(Seurat = seurat.rand, Scrna = scrna.rand, SCT = sct.rand, Ensemble = ensemble.rand))
  
  
}