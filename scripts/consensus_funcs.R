## #############################################################################
## #############################################################################
## helper functions

## identify parent node of each node in dendrogram
.pd_map <- function(hc, n) {
  ## determine parent branch node for all children nodes along dendrogram
  pd_pairs <- rbind(cbind(hc$merge[, 1], 1:(n-1)), 
                    cbind(hc$merge[, 2], 1:(n-1)))
  pd_map <- data.frame(pd_pairs[pd_pairs[, 1] > 0, ])
  names(pd_map) <- c("dtr", "prt")
  pd_map <- pd_map$prt[order(pd_map$dtr)] #the parent of each daughter
  pd_map <- c(pd_map, n) #add final node without a parent
  
  ## flip index, hclust and shc use reversed ordering
  n - rev(pd_map)
}


## determine obs indices at each node of the dendrogram
.idx_hc <- function(hc, n) {
  ## list array of cluster indices at each of the n-1 merges
  idx_hc <- array(list(), c(2*n-1, 2))
  idx_hc[1:n, 1] <- as.list(n:1)
  idx_hc[(n+1):(2*n-1), ] <- hc$merge + n + (hc$merge<0)
  
  ## complete idx_hc
  for (k in 1:(n-1)) {
    idx_hc[[n+k, 1]] <- unlist(idx_hc[idx_hc[[n+k, 1]], ])
    idx_hc[[n+k, 2]] <- unlist(idx_hc[idx_hc[[n+k, 2]], ])
  }
  
  ## flip index, hclust and shc use revered ordering
  idx_hc[(2*n-1):(n+1), ]
}



.fwer_cutoff.matrix <- function(obj, alpha, ...) {
  alpha/(nrow(obj)+1) *
    apply(obj, 1, function(x) { length(unlist(x)) })
}


.hclustFunction<-function(x,k){
   res <- PPCI::mdh(x)
   return(list(cluster = res$cluster))
}


.prepMatrix <- function(covs, ncomp = 5, ng){
  
  p <- nrow(covs[[1]])
  S <- Matrix::Matrix(0, nrow = p, ncol = p)
  stopifnot(Matrix::isSymmetric(S))
  n <- sum(ng)
  k = length(covs)
  
  
  for(i in 1:k) {
    val <- ng[i] / n
    S <- S +  (val * covs[[i]])
    gc()
  }
  
  res <- RSpectra::eigs_sym(as.matrix(S), k = ncomp)
  gc()
  
  q0 <- Matrix::Matrix(res$vectors[, 1:ncomp, drop = FALSE])
  
  return(list(A = S, v0 = q0))
  
}



#' integrate_similarity_matrices: integrate similarity matrices using a tensor product graph
#' linear combination and diffusion technique
#' 
#' Given a list of similarity matrices this function will integrate them running
#' the Shu algorithm, also can reduce noise if the input is a list consisting of
#' a single matrix.
#' 
#' @param kernellist A list of similarity matrices: those to be integrated
#' @param KNNs_p Numerical value: number of nearest neighbours for KNN graph (default=10, suggested=10-20)
#' @param diffusion_iters Numerical value: number of iterations for graph diffusion (default=4, suggested=2-6)
#' @param method Character: either TPG (see reference below) or mean (default=TPG)
#' 
#' @references Shu, Le, and Longin Jan Latecki. "Integration of single-view graphs with 
#' diffusion of tensor product graphs for multi-view spectral clustering." Asian Conference 
#' on Machine Learning. 2016.
#' 
#' @return An integrated similarity matrix
#' @export
#'
#' @examples
#' i_test <- integrate_similarity_matrices(misslfilled,method='mean')
.integrate_similarity_matrices <- function(kernellist,KNNs_p=10,diffusion_iters=4,
                                          method='TPG'){
  A <- Reduce('+', kernellist) # construct A by linear combination
  ### diffusion on KNN graph
  if (method=='TPG'){
    ## get KNN graph
    for (col in seq(1,ncol(A))){
      KNNs <- head(rev(sort(A[,col])),(KNNs_p+1)) # find the KNNs (10 default)
      tokeep <- names(KNNs)
      A[!(names(A[,col])%in%tokeep),col] <- 0
    }
    A <- A/rowSums(A) # row normalise A
    ## diffusion iterations
    Qt <- A
    im <- matrix(ncol=ncol(A),nrow=ncol(A))
    im[is.na(im)] <- 0
    diag(im) <- 1
    for (t in seq(1,diffusion_iters)){ # diffusion_iterations (4 default)
      Qt <- Matrix::crossprod(A, Matrix::crossprod(Qt, t(A)))+im
    }
    A2 <- t(Qt)
  }else if (method=='mean'){
    # if we are not doing TPG method use simple mean A
    A2 <- A/length(kernellist)
  }
  #
  return(A2)
}


.ComputeK <- function(mat, ncomp = 100){
  # result matrix
  res <- matrix(nrow=ncomp,ncol=2)
  vecs <- RSpectra::eigs_sym(as.matrix(mat), k = ncomp)$vectors
  
  # run the dip test on each component
  for (ii in seq(1,ncol(vecs))){
    r <- diptest::dip.test(vecs[,ii])
    res[ii,1] <- r$p.value
    res[ii,2] <- r$statistic
  }
}


.Reciprocal <- function(x){
  
  n = length(x)
  temp = c()
  for(i in 1: n){
    if(x[i] == 0) temp[i] = 0
    else temp[i] = 1/x[i]
  }
  return(temp)
}

.laplacian <- function(A, normalised = F){
  
  n = dim(A)[1]
  temp = apply(abs(A), 2, sum)
  D = diag(temp, nrow = n)
  
  temp1 = .Reciprocal(sqrt(temp))
  half.D = diag(temp1, nrow = n)
  if(normalised == TRUE) 	return(half.D %*% (D - A) %*% half.D)
  if(normalised == FALSE) return(D - A)
  
}


.RunConsensus <- function(obj, method = "Seurat"){
  ensemble.mat <- c()
  message("Running Significant Testing....")
  if(method == "Seurat"){
    # for Cosine
    message("Cosine....")
    eucl.eigen.labs <- .ConsensusTest(obj@Transformations$Seurat$Cosine$Eigen)
    if(length(unique(eucl.eigen.labs)) > 1){
      ensemble.mat <- cbind(ensemble.mat, eucl.eigen.labs)
    }
    # compute ensemble 
    if(!is_empty(ensemble.mat)){
      if(ncol(ensemble.mat) > 1){
        ensemble <- apply(ensemble.mat, MARGIN = 2, as.cl_partition)
        cons <- cl_consensus(as.cl_ensemble(ensemble), method = "SE", control = list(nruns = 50))
        ensemble.labs <- c(cl_class_ids(cons))
      } else if(ncol(ensemble.mat) == 1){
        ensemble.labs <- ensemble.mat[, 1]
      } 
    }
    else(
      ensemble.labs <- rep(1, nrow(obj@Transformations$Seurat$Cosine$Eigen))
    )
    # return the labels
    obj@Labels$Seurat <- ensemble.labs
  
  } else if(method == "SCTransform"){
    # for Cosine
    message("Cosine....")
    eucl.eigen.labs <- .ConsensusTest(obj@Transformations$SCTransform$Cosine$Eigen)
    if(length(unique(eucl.eigen.labs)) > 1){
      ensemble.mat <- cbind(ensemble.mat, eucl.eigen.labs)
    }
    if(!is_empty(ensemble.mat)){
      if(ncol(ensemble.mat) > 1){
        ensemble <- apply(ensemble.mat, MARGIN = 2, as.cl_partition)
        cons <- cl_consensus(as.cl_ensemble(ensemble), method = "SE", control = list(nruns = 50))
        ensemble.labs <- c(cl_class_ids(cons))
      } else if(ncol(ensemble.mat) == 1){
        ensemble.labs <- ensemble.mat[, 1]
      } 
    }
    else(
      ensemble.labs <- rep(1, nrow(obj@Transformations$SCTransform$Cosine$Eigen))
    )
    # return the labels
    obj@Labels$SCTransform <- ensemble.labs
    
  } else if(method == "Scrna"){
    # for Cosine
    message("Cosine....")

    eucl.eigen.labs <- .ConsensusTest(obj@Transformations$Scrna$Cosine$Eigen)
    if(length(unique(eucl.eigen.labs)) > 1){
      ensemble.mat <- cbind(ensemble.mat, eucl.eigen.labs)
    }
    
    if(!is_empty(ensemble.mat)){
      if(ncol(ensemble.mat) > 1){
        ensemble <- apply(ensemble.mat, MARGIN = 2, as.cl_partition)
        cons <- cl_consensus(as.cl_ensemble(ensemble), method = "SE", control = list(nruns = 50))
        ensemble.labs <- c(cl_class_ids(cons))
      } else if(ncol(ensemble.mat) == 1){
        ensemble.labs <- ensemble.mat[, 1]
      } 
    }
    else(
      ensemble.labs <- rep(1, nrow(obj@Transformations$Scrna$Cosine$Eigen))
    )
    # return the labels
    obj@Labels$Scrna <- ensemble.labs
  } else if(method == "All"){
    all.eigen.labs <- .ConsensusTest(as.matrix(obj@Transformations$All$PCA))
    
    obj@Labels$All <- all.eigen.labs
  }
  message("Done")
  return(obj)
  
}

.ConsensusTest <- function(x, alpha = 0.001, nboot = 10, nmin = 10, seed = 1994, retTree = F, method = "dc"){
  set.seed(seed)
  n <- nrow(x)
  p <- ncol(x)
  nd_type <- rep("", n-1)
  p_emp <- rep(0, n-1)
  
  
  # perform consensus clustering
  cat("\nNow Running Consensus Test\n")
  # Generate initial tree
  d <- parallelDist(scale(x), method = "euclidean")
  output <- fastcluster::hclust(d, method = "ward.D2")
  
  # do some processing
  hc_dat <- output
  idx_hc <- .idx_hc(output, n)
  cutoff <- .fwer_cutoff.matrix(idx_hc, alpha)
  pd_map <- .pd_map(output, n)
  nd_type <- rep("", n-1)
  n_min <- ncol(x)
  

  # run significance testing on each node
  for (k in 1:(n-1)) {
    ## indices for subtree
    idx_sub <- unlist(idx_hc[k, ])
    n_sub <- length(idx_sub)
    
    ## only calc p-values for branches w/ more than n_min
    if (n_sub < n_min) {
      nd_type[k] <- "n_small"
      next
    }
    
    if ((alpha < 1) && (k > 1) && (nd_type[pd_map[k]] != "sig")) {
      nd_type[k] <- "no_test"
      p_emp[k] <- 1
      next
    }
    
    # Generate initial assingments
    assignments <- kmeans(as.matrix(x[idx_sub, ]), 2, nstart = 25, iter.max = 100)$cluster
    res.pval <- UNPaC::UNPaC_Copula(x[idx_sub, ], assignments, kmeans, nsim = nboot, cov = "glasso")$pvalue_norm
    print(res.pval)
    
    if(alpha < 1){
      nd_type[k] <- ifelse(res.pval < cutoff[k],
                           "sig", "not_sig")
      p_emp[k] <- res.pval
    }
    
  }
  res <- list(in_mat = x,
              nd_type = nd_type,
              p_emp = p_emp,
              idx_hc = idx_hc,
              hc_dat = hc_dat,
              pd_map = pd_map,
              evecs = x)
  
  # get the labels
  res.labels <- .MultiCutTRee(res)
  
  if(retTree){
    return(list(Tree = res, Labels = res.labels))
  }

  return(res.labels)
  
}






.estkTW <- function(dat, ncomp = 25) {
  
  p <- ncol(dat) # number of cells
  n <- nrow(dat) # number of genes
  
  # compute Tracy-Widom bound
  x <- scale(dat)
  muTW <- (sqrt(n - 1) + sqrt(p))^2
  sigmaTW <- (sqrt(n - 1) + sqrt(p)) * (1/sqrt(n - 1) + 1/sqrt(p))^(1/3)
  sigmaHatNaive <- Matrix::crossprod(dat)
    # mixKernel::compute.kernel(x)$kernel
    # Matrix::crossprod(dat)  # x left-multiplied by its transpose
  bd <- 3.273 * sigmaTW + muTW  # 3.2730 is the p=0.001 percentile point for the Tracy-Widom distribution
  
  # compute eigenvalues and return the amount which falls above the bound
  res <- RSpectra::eigs_sym(sigmaHatNaive, k = 25, which = "LM")
  evals <- res$values
  evecs <- res$vectors
  
  # plot(evals, main = "Plot of EigenValues")
  # abline(h = bd, col = "red")
  # res <- EM_finder(evecs[, 1:sum(evals > bd )], silent = FALSE)
  # k <- findk(res, maxk = sum(evals > bd ))
  message(paste("Number of significant values:", sum(evals > bd)))
  k <- estimate_k(sigmaHatNaive, maxk = sum(evals > bd), showplots = FALSE)
  
  return(k)
  
  # return(list(numSig = sum(evals > bd ), Mat = sigmaHatNaive, Spectrum = res))
  
}


#' Calculate a distance matrix
#'
#' Distance between the cells, i.e. columns, in the input expression matrix are
#' calculated using the Euclidean, Pearson and Spearman metrics to construct
#' distance matrices.
#'
#' @param data expression matrix
#' @param method one of the distance metrics: 'spearman', 'pearson', 'euclidean'
#' @return distance matrix
#'
#' @importFrom stats cor dist
#' 
#' @useDynLib SC3
#' @importFrom Rcpp sourceCpp
#'
.calculate_distance <- function(data, method) {
  return(if (method == "spearman") {
    as.matrix(1 - cor(data, method = "spearman"))
  } else if (method == "pearson") {
    as.matrix(1 - cor(data, method = "pearson"))
  } else {
    as.matrix(parallelDist(t(data)))
  })
}

#' Distance matrix transformation
#'
#' All distance matrices are transformed using either principal component 
#' analysis (PCA) or by calculating the 
#' eigenvectors of the graph Laplacian (Spectral). 
#' The columns of the resulting matrices are then sorted in 
#' descending order by their corresponding eigenvalues.
#'
#' @param dists distance matrix
#' @param method transformation method: either 'pca' or
#' 'laplacian'
#' @return transformed distance matrix
#'
#' @importFrom stats prcomp cmdscale
#'
.transformation <- function(dists, method, ncomp = 5) {
  if (method == "pca") {
    t <- gmodels::fast.prcomp(dists, center = TRUE, scale. = TRUE)
    return(t$rotation[, 1:5])
  } else if (method == "laplacian") {
    L <- SC3::norm_laplacian(dists)
    l <- RSpectra::eigs_sym(L, k = ncomp, which = "LM")$vectors
    # sort eigenvectors by their eigenvalues
    return(l)
  }
}
