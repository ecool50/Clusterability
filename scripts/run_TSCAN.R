run_TSCAN <- function(data, trueLabels){
  
  library(TSCAN)
  data <- TSCAN::preprocess(data, takelog = TRUE, logbase = 2, cvcutoff = 0.1) #filter and normalize
  
  
  # if(clust=="set"){
  #   tsc <- TSCAN::exprmclust(data, clusternum = k_input, reduce = T)
  # }
  # 
  tsc <- TSCAN::exprmclust(data, reduce = T) #PCA

  
  clusters <- as.vector(tsc$clusterid)
  detach(package:TSCAN)
  return(list(ARI = adjustedRandIndex(trueLabels, clusters), clusters = clusters))
  
} 