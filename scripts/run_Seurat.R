run_Seurat <- function(data, ncomp = 25, trueLabels){
  
  sce <- SingleCellExperiment(
    assays = list(counts = data, logcounts = log2(data + 1)
    )
  )
  obj <- as.Seurat(
    sce,
    counts = "counts",
    data = "logcounts",
    assay = "RNA",
    project = "SingleCellExperiment"
  )
  
  obj <- Seurat::NormalizeData(object = obj, normalization.method = "LogNormalize", scale.factor = 10000)
  # VariableFeatures(obj) <- rownames(data)
  obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", verbose = F)
  obj <- Seurat::ScaleData(object = obj, features = rownames(obj@assays$RNA@data), verbose = T)
  obj <- Seurat::RunPCA(obj, features = VariableFeatures(object = obj), npcs = ncomp, verbose = F)
  obj <- Seurat::FindNeighbors(obj, reduction="pca", dims = 1:ncomp, k.param = 30)
  obj <- Seurat::FindClusters(obj, resolution = 0.8, random.seed = 1994)
  res <- Seurat::Idents(obj) 
  return(list(ARI = adjustedRandIndex(trueLabels, res), clusters = res))
}

