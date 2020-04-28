# Laplacian following Ng et al formulation
.LaplacianNg <- function(mat)
{
  D <- rowSums(mat)
  sqriD <- diag(1/sqrt(D))
  return(sqriD %*% mat %*% sqriD)
}

# Given a distance matrix, compute the gaussian similarity of its values
.distanceGaussianSimilarity <- function(distances, sigma)
{
  myfactor <- -1/(2*sigma^2)
  result   <- exp(distances^2*myfactor)
  
  # 0's are dangerous afterwards, they should be replaced by sth safer
  result[result == 0] <- 1e-16
  
  return(result)
}


# Gives the suggested sigma to spectral cluster given a distance matrix
# neighbour parameter:
# NULL -> select avg distance according to Luxburg criterium (avg distance to
# the log(number of samples)th neighbour);
# Integer > 1: avg distance to the neighbour-th closest sample
# Real number 0<neighbour<=1: avg distance to the (neighbour-th*number of samples)
# closest sample
.suggestedSigma <- function(distances, neighbour=NULL)
{
  if(is.null(neighbour))
    n <- ceiling(log(nrow(distances)))
  else #if(neighbour > 1)
    n <- neighbour
  # else
  #   n <- ceiling(neighbour*nrow(distances))
  
  dist.ord <- apply(distances, 2, sort)
  # Compute the mean removing NAs and infinite values just in case
  # If it is 0 then return 1 (or we will get errors)
  result <- mean(dist.ord[n,is.finite(dist.ord[n,])], na.rm=TRUE)
  if(result == 0.0)
    result <- 1.0
  return(result)
}


ComputeLaplacian <- function(x, sigmas = NULL, neighbours=NULL){
  nviews <- length(x)
  if(length(sigmas)     == 1) sigmas     <- rep(sigmas,     nviews)
  if(length(neighbours) == 1) neighbours <- rep(neighbours, nviews)
  
  # Placeholder to store the actual sigmas used
  mysigmas <- rep(0, nviews)
  numPoints <- nrow(as.matrix(x[[1]]))
  lapMatrix <- array(dim=c(numPoints, numPoints, nviews))
  
  for(i in 1:length(x))
  {
    if(class(x[[i]]) != "dist")
      view.dist <- as.matrix(parallelDist::parallelDist(x[[i]], method = "manhattan"))
    else
      view.dist <- as.matrix(x[[i]])
    
    if(!is.null(sigmas)) {
      mysigmas[i] <- sigmas[i]
    } else if(!is.null(neighbours)) {
      mysigmas[i] <- .suggestedSigma(view.dist, neighbours[i])
    } else {
      mysigmas[i] <- .suggestedSigma(view.dist, NULL)
    }
    view.grbf        <- .distanceGaussianSimilarity(view.dist, mysigmas[i])
    #diag(view.grbf)  <- 0
    view.lsym        <- .LaplacianNg(view.grbf)
    lapMatrix[, , i] <- view.lsym
    
    
  }
  return(lapMatrix)
  
}

# write a function to compute the top k eigenvectors
ComputeEigen <- function(lapMatrix){
  dims <- dim(lapMatrix)
  eigMatrix <- array(dim = c(dims))
  
  # compute the eigen vecotors and sort them based on eigenvalues
  for(i in 1:dims[3]){
    l <- eigen(t(lapMatrix[, , i]))
    # sort the values
    eigMatrix[, , i] <- l$vectors[, order(l$values)]
    
  }
  
  return(eigMatrix)
  
}