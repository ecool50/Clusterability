run_TSCANMod <- function(obj, trueLabels, k){
  
  library(TSCAN)
  
  # Log
  tsc <- TSCAN::exprmclust(obj@Seurat$Scaled.Data, reduce = T, k-5:k+5) #PCA
  log.clusters <- as.vector(tsc$clusterid)
  log.rand <- adjustedRandIndex(trueLabels, log.clusters)
  
  # Multinom
  tsc <- TSCAN::exprmclust(obj@Scrna$Deviances, reduce = T, k-5:k+5) #PCA
  multinom.clusters <- as.vector(tsc$clusterid)
  multinom.rand <- adjustedRandIndex(trueLabels, multinom.clusters)
  
  # NegBinom
  tsc <- TSCAN::exprmclust(obj@SCTransform$Scaled.Data, reduce = T, k-5:k+5) #PCA
  negbinom.clusters <- as.vector(tsc$clusterid)
  negbinom.rand <- adjustedRandIndex(trueLabels, negbinom.clusters)
  
  ensemble.mat <- cbind(log.clusters, multinom.clusters, negbinom.clusters)
  ensemble.labs <- diceR::LCA(ensemble.mat)
  ensemble.rand <- adjustedRandIndex(trueLabels, ensemble.labs)
  detach(package:TSCAN)
  return(list(Seurat = log.rand, Scrna = multinom.rand, SCT = negbinom.rand, Ensemble = ensemble.rand))
  
} 