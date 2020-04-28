run_CIDR <- function(data, ncomp = 5, nCluster = NULL, trueLabels){
  
  library(cidr)
  
  cid <- cidr::scDataConstructor(data, tagType = "raw") 
  cid <- cidr::determineDropoutCandidates(cid)
  cid <- cidr::wThreshold(cid)
  cid <- cidr::scDissim(cid, threads = 12, correction = FALSE)
  cid <- cidr::scPCA(cid, plotPC=FALSE)
  

  cid <- cidr::scCluster(object = cid, nCluster = NULL, nPC=ncomp) 
  
  detach(package:cidr)
  clusters <- cid@clusters
  return(list(ARI = adjustedRandIndex(trueLabels, clusters), clusters = clusters))
}
