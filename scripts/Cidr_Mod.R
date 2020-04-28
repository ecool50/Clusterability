run_AscendMod <- function(obj, truelabels, ncomp.seurat, ncomp.scrna, ncomp.sct){
  # Seurat
  em <- ascend::EMSet(list(counts = obj@Seurat$Scaled.Data, logcounts = obj@Seurat$Scaled.Data, normcounts = obj@Seurat$Scaled.Data))
  
  em <- ascend::runPCA(em, ngenes = 500, scaling = FALSE)
  
  em1 <- ascend::runCORE(em, dims = ncomp.seurat)
  
  seurat.clusters <- ascend::clusterAnalysis(em1)$clusters
  seurat.rand <- adjustedRandIndex(trueLabels, seurat.clusters)
    
  # Scrna
  em <- ascend::EMSet(list(counts = obj@Scrna$Deviances, logcounts = obj@Scrna$Deviances, normcounts = obj@Scrna$Deviances))
  
  em <- ascend::runPCA(em, ngenes = 500, scaling = FALSE)
  
  em1 <- ascend::runCORE(em, dims = ncomp.seurat)
  
  scrna.clusters <- ascend::clusterAnalysis(em1)$clusters
  scrna.rand <- adjustedRandIndex(trueLabels, scrna.clusters)
  
  # SCT
  em <- ascend::EMSet(list(counts = obj@SCTransform$Scaled.Data, logcounts = obj@SCTransform$Scaled.Data, normcounts = obj@SCTransform$Scaled.Data))
  
  em <- ascend::runPCA(em, ngenes = 500, scaling = FALSE)
  
  em1 <- ascend::runCORE(em, dims = ncomp.sct)
  
  sct.clusters <- ascend::clusterAnalysis(em1)$clusters
  sct.rand <- adjustedRandIndex(trueLabels, sct.clusters)
  
  
  return(list(Seurat = seurat.rand, Scrna = scrna.rand, SCT = sct.rand))
}