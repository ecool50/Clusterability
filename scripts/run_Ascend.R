# Use capture.output({ }) to not show package messages

run_ascend <- function(data, ncomp = 5, trueLabels, ngenes = 1500){
  
  library(ascend)
  
  em <- ascend::EMSet(list(counts = data))
  

  
  # em <- ascend::filterLowAbundanceGenes(em) #filters
  em <- ascend::normaliseByRLE(em) #normalizes 
  
  em <- ascend::runPCA(em, ngenes = ngenes)
  
  em1 <- ascend::runCORE(em, dims = ncomp)
  
  
  # else if(use_dims_input=="internal"){
  #   em1 <- ascend::runCORE(em, nres=40)
  # }
  
  clusters <- ascend::clusterAnalysis(em1)$clusters
  missing_ids=which(colnames(em) %in% setdiff(colnames(em), colnames(em1)))
  
  detach(package:ascend)
  
  if(length(missing_ids)==0){
    return(list(ARI = adjustedRandIndex(trueLabels, clusters), clusters = clusters))
  } else{
    return(list(clusters=clusters, missing_ids=missing_ids))
  }
  
  
}






