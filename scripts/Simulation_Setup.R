Setup_1 <- function(numCells = 5000, numGenes = 1000, numClusters = 4, size = "balanced"){
  suppressPackageStartupMessages(library(splatter))
  # balanced cluster sizes
  if(size == "balanced"){
    group.prob <- rep(1/numClusters, numClusters)
  } else if(size == "unbalanced"){
    x <- runif(numClusters)
    group.prob=x/sum(x)
  }
  
  # set routine parameters
  de.prob=0.5
  dropout.type ="experiment"
  dropout.mid=0
  print(group.prob)
  # create a new splat object
  params <- newSplatParams()
  # Generate the data
  params <- setParams(params, batchCells=numCells, nGenes=numGenes, 
                      group.prob = unlist(group.prob), de.prob=de.prob, dropout.type=dropout.type, seed=1994)
  sce <- splatSimulateGroups(params, verbose = TRUE)
  
}

Setup_2 <- function(numCells = 5000, numGenes = 1000, numClusters = 4, separability = 0){
  group.prob <- rep(1/numClusters, numClusters)
  de.prob <- separability
  dropout.type <-"experiment"
  dropout.mid <- 0
  # de.facScale <- rep(1,numClusters)*separability
  
  # create a new splat object
  params <- newSplatParams()
  # Generate the data
  params <- setParams(params, batchCells=numCells, nGenes=numGenes, 
                      group.prob = unlist(group.prob), dropout.type=dropout.type, de.prob = de.prob, seed=1994)
  sce <- splatSimulateGroups(params, verbose = TRUE)
  
  return(sce)
  
}

# setup 3
Setup_3 <- function(numCells = 5000, numGenes = 1000, numClusters = 4, dropout.mid = 0){
  de.prob=0.5
  dropout.type ="experiment"
  group.prob <- rep(1/numClusters, numClusters)
  # create a new splat object
  params <- newSplatParams()
  # Generate the data
  params <- setParams(params, batchCells=numCells, nGenes=numGenes, 
                      group.prob = unlist(group.prob), de.prob=de.prob, dropout.type=dropout.type, dropout.mid = dropout.mid, seed=1994)
  sce <- splatSimulateGroups(params, verbose = TRUE)
  return(sce)
}

# this function takes in a counts matrix and returns a modality object with tests slot filled
RunModalityTest <- function(Counts){
  message("Preprocessing the counts")
  obj <- CreateModalityObject(Counts, sparse = FALSE)
  # add the seurat slot
  message("Log Normalisation")
  obj <- PreprocessObject(obj,  nFeatures = 500)
  # add the sctransform slot
  message("Negative Binomial Normalisation")
  obj <- PreprocessObject(obj,  method = "SCTransform", nFeatures = 500)
  # add the scrna slot
  message("Multinomial Normalisation")
  obj <- PreprocessObject(obj, method = "Scrna", nFeatures = 500)
  
  # run modality testing
  message("Running the clusterability test")
  obj <- TestModality(obj)
  obj <- TestModality(obj, type = "Scrna")
  obj <-TestModality(obj, type = "SCTransform")
  
  return(unlist(obj@Tests))
  
  
}

RunSignificanceTest <- function(Counts){
  obj <- CreateModalityObject(Counts, sparse = FALSE)
  # add the seurat slot
  message("Log Normalisation")
  obj <- PreprocessObject(obj,  nFeatures = 500)
  # add the sctransform slot
  message("Negative Binomial Normalisation")
  obj <- PreprocessObject(obj,  method = "SCTransform", nFeatures = 500)
  # add the scrna slot
  message("Multinomial Normalisation")
  obj <- PreprocessObject(obj, method = "Scrna", nFeatures = 500)
  
  message("Running the significance test for clusterability")
  
  obj <- ComputeTransformations(obj)
  obj <- ComputeTransformations(obj, method = "Scrna")
  obj <- ComputeTransformations(obj, method = "SCTransform")
  obj <- ComputeTransformations(obj, method = "All")
  
  obj <- .RunConsensus(obj)
  obj <- .RunConsensus(obj, method = "Scrna")
  obj <- .RunConsensus(obj, method = "SCTransform")
  obj <- .RunConsensus(obj, method = "All")
  
  
  return(obj)
}