run_pcaReduce <- function(data, q = 5){
  
  library(pcaReduce)
  
  
  pcar <- pcaReduce::PCAreduce(t(data), nbt = 1, q = q  , method = 'M')
  clusters <- as.vector(pcar)[[1]][,1]
  # detach(package:pcaReduce)
  return(clusters)
}