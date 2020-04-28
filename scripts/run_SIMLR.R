run_SIMLR <- function(data ,dims_input = 25, clust = "estimate", k, trueLabels){
  
  library(SIMLR)
  
  
  # if(clust=="estimate"){
  #   est = SIMLR::SIMLR_Estimate_Number_of_Clusters(data, NUMC = k-2:k+2) #2:5 is default
  #   k_input = (k-2:k+2)[which.min(est$K1)]
  # }
  
  # if(preproc=="none" && use_dims_input=="TRUE"){
  #   siml <- SIMLR::SIMLR(X = data, c = k_input, normalize = FALSE, no.dim=dims_input)
  # }
  
  # else if(preproc=="none" && use_dims_input=="internal"){
  #   siml <- SIMLR::SIMLR(X = data, c = k_input, normalize = FALSE, no.dim=NA)
  # }
  
  
  siml <- SIMLR::SIMLR_Large_Scale(X = data, c = k, normalize = TRUE)
  
  
  # else if(preproc=="method_specific" && use_dims_input=="internal"){
  #   siml <- SIMLR::SIMLR(X = data, c = k_input, normalize = TRUE, no.dim=NA)
  # }
  
  clusters <- siml$y$cluster
  detach(package:SIMLR)
  return(list(ARI = adjustedRandIndex(trueLabels, clusters), clusters = clusters))
}