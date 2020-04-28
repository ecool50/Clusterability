run_SeuratMod <- function(obj, ncomp.seurat, ncomp.scrna, ncomp.sct, trueLabels){
  # Seurat
  seurat.obj <- obj@Seurat$obj
  seurat.obj <- Seurat::FindNeighbors(seurat.obj, dims = 1:5)
  seurat.obj <- Seurat::FindClusters(seurat.obj, resolution = 0.8, random.seed = 1994)
  res.seurat <- Seurat::Idents(seurat.obj) 
  seurat.rand <- adjustedRandIndex(trueLabels, res.seurat)
  
  # Scrna
  sce <- SingleCellExperiment(
    assays = list(counts = obj@Scrna$Deviances, logcounts = obj@Scrna$Deviances
    )
  )
  obj.scrna <- as.Seurat(
    sce,
    counts = "counts",
    data = "logcounts",
    assay = "RNA",
    project = "SingleCellExperiment"
  )
  
  obj.scrna@assays$RNA@scale.data <- as.matrix(obj@Scrna$Deviances)
  VariableFeatures(obj.scrna) <- rownames(obj@Scrna$Deviances)
  obj.scrna <- Seurat::RunPCA(obj.scrna, features = VariableFeatures(object = obj.scrna), npcs = ncomp.scrna, verbose = F)
  obj.scrna <- Seurat::FindNeighbors(obj.scrna, dims = 1:ncomp.scrna)
  obj.scrna <- Seurat::FindClusters(obj.scrna, resolution = 0.8, random.seed = 1994)
  res.scrna <- Seurat::Idents(obj.scrna) 
  scrna.rand <- adjustedRandIndex(trueLabels, res.scrna)
  
  
  # SCTransform
  sce <- SingleCellExperiment(
    assays = list(counts = obj@SCTransform$Scaled.Data, logcounts = obj@SCTransform$Scaled.Data
    )
  )
  obj.sct <- as.Seurat(
    sce,
    counts = "counts",
    data = "logcounts",
    assay = "RNA",
    project = "SingleCellExperiment"
  )
  
  obj.sct@assays$RNA@scale.data <- as.matrix(obj@SCTransform$Scaled.Data)
  VariableFeatures(obj.sct) <- rownames(obj@SCTransform$Scaled.Data)
  obj.sct <- Seurat::RunPCA(obj.sct, features = VariableFeatures(object = obj.sct), npcs = ncomp.sct, verbose = F)
  obj.sct <- Seurat::FindNeighbors(obj.sct, dims = 1:ncomp.sct)
  obj.sct <- Seurat::FindClusters(obj.sct, resolution = 0.8, random.seed = 1994)
  res.sct <- Seurat::Idents(obj.sct) 
  sct.rand <- adjustedRandIndex(trueLabels, res.sct)
  
  
    # compute ensemble
  ensemble.mat <- cbind(res.seurat, res.scrna, res.sct)
  ensemble.labs <- diceR::LCA(ensemble.mat)
  ensemble.rand <- adjustedRandIndex(trueLabels, ensemble.labs)
  return(list(Seurat = seurat.rand, Scrna = scrna.rand, SCT = sct.rand, Ensemble = ensemble.rand))
}

