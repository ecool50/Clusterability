
# This script is an attempt to implement the clusterability pipeline
# 
# Elijah Willie
# July 08 2019


# Create a function to preprocess the data using either the standard Seurat pipeline the Multinomial modelling of the 
# raw counts as done in https://www.biorxiv.org/content/biorxiv/early/2019/03/11/574574.full.pdf

source("Helper_Functions.R")
options(warn=-1)

# write a function to create a modality object
CreateModalityObject <-function(counts, sparse = TRUE){
  # create a generator function
  setClass("Modality", slots = c("Data", "Seurat", "SCTransform", "Randomly", "Scrna", "Dists", "Transformations", "Tests", "Plots", "Labels"))
  
  # generate the modality object
  obj <- new("Modality")
  obj@Data$counts <- Matrix(counts, sparse = sparse)
  
  # return the obj
  return(obj)
  
}

PreprocessObject <- function(object, method = "Seurat", nPCs = 100, nFeatures = 500, prop = 0.10){
  # dat is a counts matrix where the cells are columns and genes are rows
  
  
  if(!(method %in% c("Seurat", "Scrna", "SCTransform", "Randomly"))){
    stop(paste(type, "is not a valid option"))
  }
  
  # if we are preprocessing using Seurat
  if(method == "Seurat"){
    
    message("\nRunning Log Preprocessing\n")
    # create the seurat object
    obj <- Seurat::CreateSeuratObject(counts = object@Data$counts, project = "Proj", min.cells = 5, min.features = 200)
    
    # compute % mitohondrial genes
    obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^mt-")
    # subset the data
    obj <- subset(obj, subset = nFeature_RNA > 0 & percent.mt < 5)
    cleaned_data <- as.matrix(obj@assays$RNA@data)
    
    # normalize and scale the data
    obj <- Seurat::NormalizeData(object = obj, normalization.method = "LogNormalize", verbose = F)
    
    # compute highly variable genes
    obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = nFeatures, verbose = F)
    
    # remove unwanted sources of variation
    obj <- Seurat::ScaleData(object = obj, features = VariableFeatures(obj), verbose = T)
    
    # extract the scaled data
    scaled_data <- obj@assays$RNA@scale.data
    
    message("Estimating Significant PCs")
    # nPCs <- round(intrinsicDimension::maxLikGlobalDimEst(t(scaled_data), unbiased = TRUE, k = 15)$dim.est)
    
    if (ncol(scaled_data) > 5000){
      inds <- sample(ncol(scaled_data), 5000)
      x_train <- scaled_data[, inds]
      
      ncomp <-  min(10,.SigPCs(x_train, B = 50)$r)
      pca <- prcomp(t(x_train))
      pca.pred <- predict(pca, t(scaled_data))
      sig_pcs_data <- pca.pred[, 1:max(2,ncomp)]
    } else{
      sig_pcs_data <-  .SigPCs(scaled_data, B = 50)$sigVecs
    }
    
    # nPCs <- max(5, nPCs)
    
    # run PCA on the highly variable genes
    # sig_pcs_data <- irlba::prcomp_irlba(scaled_data, n = nPCs, retx = T, fastpath = T)$rotation
    
    
    # Run Tsne
    message("\nRunning Tsne")
    # obj <- Seurat::RunTSNE(obj, dims = 1:nPCs, check_duplicates=FALSE)
    tsne.res <- as.matrix(Rtsne::Rtsne(t(scaled_data), perplexity=30, verbose=F, max_iter = 500, pca = T,
                                       normalize = T, check_duplicates=FALSE, num_threads = 0)$Y)
    rownames(tsne.res) <- colnames(scaled_data)
    # 
    # # Run UMAP
    # message("\nRunning UMAP")
    # # obj <- Seurat::RunUMAP(obj, dims = 1:nPCs, n.components = nPCs)
    # umap.res <- uwot::umap(t(scaled_data), n_neighbors = 30, metric = "cosine", min_dist = 0.3, n_threads = 12)

    # message("Computing Principal Curve")
    # fit <- pcurve(sig_pcs_data, robust = FALSE, plot.pca = F, plot.true = F, plot.init = F,
    #               plot.segs = F, plot.resp = F, plot.cov = F, start = sig_pcs_data[,1])
    
    # update the object and return it
    object@Data$counts <- NULL
    object@Data$cleaned <- Matrix(cleaned_data, sparse = TRUE)
    object@Seurat$Scaled.Data <- Matrix(scaled_data, sparse = TRUE)
    object@Seurat$Sig.PCs.Data <- sig_pcs_data
    object@Seurat$Tsne.Data <- tsne.res
    # object@Seurat$Umap.Data <- umap.res
    # # object@Seurat$Autoencoder <- dat.encoded
    # # object@Seurat$Princurve <- fit
    # object@Seurat$obj <- obj
    
   # object@Seurat$Phate.Data <- phate.res
    return(object)
  }
  
  
  # if we are preprocessing using Seurat
  if(method == "SCTransform"){
    
    message("\nRunning NegBinom Preprocessing\n")
    data <- object@Data$cleaned
    future::plan(strategy = 'multicore')
    if(ncol(data) > 5000){
      ncells <- 5000
    } else{
      ncells <- ncol(data)
    }
    vst_out <- sctransform::vst(round(data), latent_var = c('log_umi'), return_gene_attr = TRUE, return_cell_attr = TRUE, 
                                show_progress = TRUE, residual_type = "pearson", n_cells = ncells)
    
    # get the most deviant genes
    top_genes <- rownames(round(vst_out$gene_attr[order(-vst_out$gene_attr$residual_variance), ], 2)[1:nFeatures, ])
    
    # subset matrix
    scaled_data <- vst_out$y[top_genes, ]
    
    message("Computing Significant PCs")
    # nPCs <- round(intrinsicDimension::maxLikGlobalDimEst(t(scaled_data), unbiased = TRUE, k = 15)$dim.est)
  
    # nPCs <- max(5, nPCs)
    # sig_pcs_data <- irlba::prcomp_irlba(scaled_data, n = nPCs, retx = T, fastpath = T)$rotation
    if (ncol(scaled_data) > 5000){
      inds <- sample(ncol(scaled_data), 5000)
      x_train <- scaled_data[, inds]
      
      ncomp <-  min(10,.SigPCs(x_train, B = 50)$r)
      sig_pcs_data <- irlba::prcomp_irlba(scaled_data, n = max(3,ncomp))$rotation
      
    } else{
      sig_pcs_data <-  .SigPCs(scaled_data, B = 50)$sigVecs
    }
    # Run Tsne
    message("\nRunning Tsne")
    tsne.res <- as.data.frame(Rtsne::Rtsne(t(scaled_data), perplexity=30, verbose=F, max_iter = 500, pca = T,
                                       normalize = T, check_duplicates=FALSE, num_threads = 0)$Y)
    rownames(tsne.res) <- colnames(scaled_data)
    # 
    # # Run UMAP
    # message("\nRunning UMAP")
    # # obj <- Seurat::RunUMAP(obj, dims = 1:nPCs, n.components = nPCs)
    # umap.res <- uwot::umap(t(scaled_data), n_neighbors = 30, metric = "cosine", min_dist = 0.3, n_threads = 12)
    
    # message("Computing Principal Curve")
    # fit <- pcurve(sig_pcs_data, robust = FALSE, plot.pca = F, plot.true = F, plot.init = F,
    #               plot.segs = F, plot.resp = F, plot.cov = F, start = sig_pcs_data[, 1])
  
    # update the object and return it
    object@SCTransform$Scaled.Data <- Matrix(scaled_data, sparse = TRUE)
    object@SCTransform$Sig.PCs.Data <- sig_pcs_data
    object@SCTransform$Tsne.Data <- tsne.res
    # object@SCTransform$Umap.Data <- umap.res
    # object@SCTransform$Autoencoder <- dat.encoded
    # object@SCTransform$Princurve <- fit
    # object@Seurat$Phate.Data <- phate.res
    return(object)
  }
  
  
  # if we are preprocessing using multinomial modelling
  if(method == "Scrna"){
    
    message("\nRunning Multinomial Preprocessing\n")
    dat <- object@Data$cleaned
    
    # Get the most informative genes
    message("Computing Informative Genes")
    d  <- scry::devianceFeatureSelection(as.matrix(dat))
    d <- data.frame(deviance = d)
    rownames(d) <- rownames(dat)
    x <-order(d$deviance, decreasing = TRUE)[1:nFeatures]
    top_genes <- rownames(d)[x]
    
    gc()
    
    # extract the counts containing the top n most informative genes
    counts_filtered <- dat[top_genes, ]
    gc()
    
    # compute the poisson deviances
    message("Computing Null Residuals")
    obj_resids <- scry::nullResiduals(counts_filtered, fam = "poisson")
    gc()
    
    message("Computing Significant PCs")
    # nPCs <- round(intrinsicDimension::maxLikGlobalDimEst(t(obj_resids), unbiased = TRUE, k = 15)$dim.est)
    # nPCs <- max(5, nPCs)
    
    # obj_PCA <- irlba::prcomp_irlba(obj_resids, n = nPCs, retx = T, fastpath = T)$rotation
    if (ncol(obj_resids) > 5000){
      inds <- sample(ncol(obj_resids), 5000)
      x_train <- obj_resids[, inds]
      
      ncomp <- min(10,.SigPCs(x_train, B = 50)$r)
      sig_pcs_data <- irlba::prcomp_irlba(obj_resids, n = max(3,ncomp))$rotation
    } else{
      sig_pcs_data <-  .SigPCs(obj_resids, B = 50)$sigVecs
    }
    # obj_PCA <- scry::GLMPCA(as.matrix(counts_filtered), L = nPCs, fam = "mult", verbose = TRUE, ctl = list(maxIter = 100, eps = 1e-04))
    
    # Run Tsne
    message("\nRunning Tsne")
    tsne.res <- as.matrix(Rtsne::Rtsne(t(obj_resids), perplexity=30, verbose=F, max_iter = 500, pca = T,
                                       normalize = T, check_duplicates=FALSE, num_threads = 0)$Y)
    rownames(tsne.res) <- colnames(obj_resids)
    # 
    # # Run UMAP
    # message("\nRunning UMAP")
    # umap.res <- uwot::umap(t(obj_resids), n_neighbors = 30, metric = "cosine", min_dist = 0.3, n_threads = 12)
    
    # message("Computing Principal Curve")
    # fit <- pcurve(sig_pcs_data, robust = FALSE, plot.pca = F, plot.true = F, plot.init = F,
    #               plot.segs = F, plot.resp = F, plot.cov = F, start = sig_pcs_data[, 1])
    
    # fit <- princurve::principal_curve(sig_pcs_data)
    
    #update the object and return it
    object@Scrna$Deviances <- Matrix(obj_resids, sparse = TRUE)
    object@Scrna$Sig.PCs.Data <- sig_pcs_data
    object@Scrna$Tsne.Data <- tsne.res
    # object@Scrna$Umap.Data <- umap.res
    # object@Scrna$Autoencoder <- dat.encoded
    # object@Scrna$Princurve <- fit
    
    return(object)
  }
}

# Create a function to test multimodality

TestModality <- function(obj, type = "Seurat"){
  
  if(!(type %in% c("Scrna", "Seurat", "SCTransform", "Randomly"))){
    stop(paste(type, "is not a valid option"))
  }
  
  # if we are dealing with PC data
  if(type == "Seurat"){
    # Dip test
    message('Running Dip Test')
    res.pca <- .RunDip(obj@Seurat$Sig.PCs.Data, "PCA")
    # res.prin <- .RunDip(obj@Seurat$Princurve$lambda, "Princurve")
    
    # Silverman
    # message('Running Silverman Test')
    # res.pca.silver <- .RunSilver(obj@Seurat$Sig.PCs.Data[, 1], "PCA")
    # res.prin.silver <- .RunSilver(obj@Seurat$Princurve$lambda, "Princurve")
    # compile the results and update the object
    res <- list("PCA" = res.pca)
    # res.sliver <- list("PCA" = res.pca.silver, "Princurve" = res.prin.silver)
    obj@Tests$Seurat <- list(Dip = res)
   
    return(obj)
  }

  if(type == "Scrna"){
    # Diptest
    message('Running Dip Test')
    res.pca <- .RunDip(obj@Scrna$Sig.PCs.Data, "PCA")
    # res.prin <- .RunDip(obj@Scrna$Princurve$lambda, "Princurve")
    
    # Silverman
    # message('Running Silverman Test')
    # res.pca.silver <- .RunSilver(obj@Scrna$Sig.PCs.Data[, 1], "PCA")
    # res.prin.silver <- .RunSilver(obj@Scrna$Princurve$lambda, "Princurve")
    # compile the results and update the object
    res <- list("PCA" = res.pca)
    # res.sliver <- list("PCA" = res.pca.silver, "Princurve" = res.prin.silver)
    obj@Tests$Scrna <- list(Dip = res)
    
    return(obj)
  }
  
  if(type == "SCTransform"){
    # Diptest
    message('Running Dip Test')
    res.pca <- .RunDip(obj@SCTransform$Sig.PCs.Data, "PCA")
    # res.prin <- .RunDip(obj@SCTransform$Princurve$lambda, "Princurve")
    
    # # Silverman
    # message('Running Silverman Test')
    # res.pca.silver <- .RunSilver(obj@SCTransform$Sig.PCs.Data[, 1], "PCA")
    # res.prin.silver <- .RunSilver(obj@SCTransform$Princurve$lambda, "Princurve")
    # compile the results and update the object
    res <- list("PCA" = res.pca)
    # res.sliver <- list("PCA" = res.pca.silver, "Princurve" = res.prin.silver)
    obj@Tests$SCTransform <- list(Dip = res)
    
    return(obj)
  }
  
}


# create a function to run dip
.RunDip <- function(dat, type = "PC"){
  if(!(type %in% c( "PCA", "Multiview", "Princurve"))){
    stop(paste(type, "is not a valid option"))
  }
    
    if(type == "Princurve"){
      message("\nRunning tests on Principal Curve data\n")
      pval <- dip.test(dat)$p.value
      return(pval)
    } 
    
    # if we are using the projected data
    if(type == "MDH"){
      message("\nRunning tests on MDH data\n")
      pval <- diptest::dip.test(dat)$p.value
      # return the results
      return(pval)
      
    }
  
   # if we are using Tsne
   if(type == "PCA"){
     message("\nRunning tests on PCA data\n")
     dcos <- as.dist(coop::cosine(t(dat)))
     pval.cos <- diptest::dip.test(dcos)$p.value
     
     gc()
     return(pval.cos)
     
   }
  
  if(type == "Multiview"){
    cat("\nRunning tests on Multiview data\n")
    pval <- dip.test(dat)$p.value
    return(pval)
  }
}


# create a function to run dip
.RunSilver <- function(dat, type = "PCA"){
  if(!(type %in% c( "PCA", "Multiview", "Princurve"))){
    stop(paste(type, "is not a valid option"))
  }
  
  if(type == "Princurve"){
    message("\nRunning tests on Principal Curve data\n")
    pval <- multimode::modetest(dat, method = "SI", B = 100)$p.value
    return(pval)
  } 
  
  # if we are using the projected data
  if(type == "MDH"){
    message("\nRunning tests on MDH data\n")
    pval <- multimode::modetest(dat, method = "SI", B = 100)$p.value
    # return the results
    return(pval)
    
  }
  
  # if we are using Tsne
  if(type == "PCA"){
    message("\nRunning tests on PCA data\n")
    pval <- multimode::modetest(dat, method = "SI", B = 100)$p.value
    gc()
    return(pval)
    
  }
  
  if(type == "Multiview"){
    cat("\nRunning tests on Multiview data\n")
    pval <- multimode::modetest(dat, method = "SI", B = 100)$p.value
    return(pval)
  }
}


.SigPCs <- function(dat,
                          B = 100, threshold = 0.05,
                          verbose = T, seed = NULL, ndf = 25) {
  if (!is.null(seed))
    set.seed(seed)
  n <- ncol(dat)
  m <- nrow(dat)
  
  uu <- RSpectra::svds(dat, k = ndf, nu = 0)
  dstat <- uu$d[1:ndf]^2/sum(uu$d[1:ndf]^2)
  dstat0 <- matrix(0, nrow = B,
                   ncol = ndf)
  if (verbose == TRUE)
    message("\nEstimating a number of significant principal component: ")
  for (i in 1:B) {
    if (verbose == TRUE)
      cat(paste(i, " "))
    dat0 <- t(apply(dat, 1,
                    sample, replace = FALSE))
    uu0 <-  RSpectra::svds(dat0, k = ndf, nu = 0)
    dstat0[i, ] <- uu0$d[1:ndf]^2/sum(uu0$d[1:ndf]^2)
  }
  p <- rep(1, ndf)
  # print(dstat0)
  for (i in 1:ndf) {
    p[i] <- mean(dstat0[, i] >=
                   dstat[i])
  }
  for (i in 2:ndf) {
    p[i] <- max(p[(i - 1)],
                p[i])
  }
  r <- sum(p < threshold)
  sig_indices <- which(p < 0.05)
  if(length(sig_indices) < 3){
    sigVecs <- uu$v[, 1:3]
  } else{
    sigVecs <- uu$v[, sig_indices]
  }
  
  return(list(r = r, p = p, sigVecs = sigVecs))
}

# function to compute distances
ComputeDistances <- function(obj, method = "Seurat"){
  message("Calculating Distances...")
  if(method == "Seurat"){
    dcor <- 1 - coop::pcor(obj@Seurat$Scaled.Data)
    deucl <- as.matrix(parallelDist::parallelDist(t(obj@Seurat$Scaled.Data)))
    dspear <- 1 - WGCNA::cor(obj@Seurat$Scaled.Data, method = "spearman", nThreads = 12)
    obj@Dists$Seurat <- list(Pearson = dcor, Euclidean = deucl, Spearman = dspear)
    
  } else if(method == "SCTransform"){
    dcor <- 1 - coop::pcor(obj@SCTransform$Scaled.Data)
    deucl <- as.matrix(parallelDist::parallelDist((t(obj@SCTransform$Scaled.Data))))
    dspear <- 1 - WGCNA::cor(obj@SCTransform$Scaled.Data, method = "spearman", nThreads = 12)
    obj@Dists$SCTransform <- list(Pearson = dcor, Euclidean = deucl, Spearman = dspear)
    
  } else if(method == "Scrna"){
    dcor <- 1 - coop::pcor(obj@Scrna$Deviances)
    deucl <- as.matrix(parallelDist::parallelDist(t(obj@Scrna$Deviances)))
    dspear <- 1 - WGCNA::cor(obj@Scrna$Deviances, method = "spearman", nThreads = 12)
    obj@Dists$Scrna <- list(Pearson = dcor, Euclidean = deucl, Spearman = dspear)
  }
  message("Done")
 return(obj)     
}

# function to compute transformations
ComputeTransformations <- function(obj, method = "Seurat"){
  message("Transforming Distances...")
  if(method == "Seurat"){
    # cosine
    x <- as.matrix(obj@Seurat$Scaled.Data)
    
    dcos <- (1 - coop::cosine(x))/2


    # dcos <- sqrt(2*dcos)
    gc()
    # ncomp <- dim(modes::amps(parallelDist(as.matrix(obj@Seurat$Princurve$lambda)))$Peaks)[1]
    # ncomp <- dim(modes::amps(dcos)$Peaks)[1]
    # ncomp <- round(intrinsicDimension::maxLikGlobalDimEst(x, unbiased = TRUE, k = 15)$dim.est)
    ncomp <- dim(obj@Seurat$Sig.PCs.Data)[2]
    ncomp <- min(ncomp,10)
    tpca <- irlba::prcomp_irlba(dcos, n = max(2, ncomp), center = T, scale. = T)$rotation
    # tpca <- RSpectra::eigs_sym(dcos, k = ncomp)$vectors
    rownames(tpca) <- colnames(obj@Seurat$Scaled.Data)
    obj@Transformations$Seurat$Cosine <- list(Eigen = tpca)
    
    } else if(method == "SCTransform"){
      
      x <- as.matrix(obj@SCTransform$Scaled.Data)
      dcos <- (1 - coop::cosine(x))/2

      # dcos <- sqrt(2*dcos)
      gc()
      # ncomp <- dim(modes::amps(parallelDist(as.matrix(obj@Seurat$Princurve$lambda)))$Peaks)[1]
      # ncomp <- dim(modes::amps(dcos)$Peaks)[1]
      # ncomp <- round(intrinsicDimension::maxLikGlobalDimEst(x, unbiased = TRUE, k = 15)$dim.est)
      ncomp <- dim(obj@Seurat$Sig.PCs.Data)[2]
      ncomp <- min(ncomp,10)
      # tpca <- RSpectra::eigs_sym(dcos, k = ncomp)$vectors
      tpca <- irlba::prcomp_irlba(dcos, n = max(2, ncomp), center = T, scale. = T)$rotation
      rownames(tpca) <- colnames(obj@SCTransform$Scaled.Data)
      obj@Transformations$SCTransform$Cosine <- list(Eigen = tpca)
      gc()
  
    } else if(method == "Scrna"){
      
      # cosine
      x <- as.matrix(obj@Scrna$Deviances)
      dcos <- (1 - coop::cosine(x))/2

      # dcos <- sqrt(2*dcos)
      gc()
      # ncomp <- dim(modes::amps(parallelDist(as.matrix(obj@Seurat$Princurve$lambda)))$Peaks)[1]
      # ncomp <- dim(modes::amps(dcos)$Peaks)[1]
      # ncomp <- round(intrinsicDimension::maxLikGlobalDimEst(x, unbiased = TRUE, k = 15)$dim.est)
      ncomp <- dim(obj@Seurat$Sig.PCs.Data)[2]
      ncomp <- min(ncomp,10)
      tpca <- irlba::prcomp_irlba(dcos, n = max(2, ncomp), center = T, scale. = T)$rotation
      rownames(tpca) <- colnames(obj@Scrna$Deviances)
      # tpca <- RSpectra::eigs_sym(dcos, k = ncomp)$vectors
      obj@Transformations$Scrna$Cosine <- list(Eigen = tpca)
      gc()
    
    } else if(method == "All"){
      
      d1 <- coop::covar(as.matrix(obj@Seurat$Scaled.Data))
      gc()
      d2 <- coop::covar(as.matrix(obj@Scrna$Deviances))
      gc()
      d3 <- coop::covar(as.matrix(obj@SCTransform$Scaled.Data))
      gc()
      
      ng <- c(nrow(d1), nrow(d2), nrow(d3))
      mat <- list(d1,d2, d3)
      A <- Reduce('+', mat)*(1/2)
      # ncomp <- dim(modes::amps(A)$Peaks)[1]
      # ncomp <- round(intrinsicDimension::maxLikGlobalDimEst(A, unbiased = TRUE, k = 15)$dim.est)
      ncomp <- min(dim(obj@Seurat$Sig.PCs.Data)[2], dim(obj@Scrna$Sig.PCs.Data)[2], dim(obj@SCTransform$Sig.PCs.Data)[2])
      tpca <- irlba::prcomp_irlba(A, n = max(2, ncomp), center = T, scale. = T)$rotation
      rownames(tpca) <- colnames(obj@Scrna$Deviances)
      # tpca <- RSpectra::eigs_sym(A, k = ncomp)$vectors
      obj@Transformations$All$PCA <- tpca
      
      gc()
    } 
  
  message("Done...")
  return(obj)
 
}

# write a function to run teh various methods
RunMethods <- function(obj){
  # compute labels
  message("\nRunning Significance Clustering On Log Model")
  obj <- .RunConsensus(obj)
  message("\nRunning Significance Clustering On Multinomial Model")
  obj <- .RunConsensus(obj, method = "Scrna")
  message("\nRunning Significance Clustering On Negative Binomial Model")
  obj <- .RunConsensus(obj, method = "SCTransform")
  # message("\nRunning Significance Clustering on Randomly Model")
  # obj <- .RunConsensus(obj, method = 'Randomly')
  message("\nRunning Significance Clustering On Multiview Model")
  obj <- .RunConsensus(obj, method = "All")
  
  # Run theother methods
  gene.names <- union(union(rownames(obj@Scrna$Deviances), rownames(obj@Seurat$Scaled.Data)), rownames(obj@SCTransform$Scaled.Data))
  data <- as.matrix(obj@Data$cleaned)[gene.names, ]
  message("\nRunning CIDR")
  res.cidr <- run_CIDR(data, trueLabels = trueLabels)
  message("\nNow Running Monocle")
  res.monnocle <- run_monocle3(data, trueLabels = trueLabels)
  # message("\nRunning SIMLR")
  # res.simlr <- run_SIMLR(data, trueLabels = trueLabels, k = length(unique(trueLabels)))
  message("\nRunning sscClust")
  res.sscClust <- run_sscClust(obj@Seurat$Scaled.Data, trueLabels = trueLabels, k = length(unique(trueLabels)))
  gc()
  message("\nRunning TSCAN")
  res.tscan <- run_TSCAN(data, trueLabels)
  gc()
  message("\nRunning SC3")
  res.sc3 <- run_SC3(data, trueLabels = trueLabels)
  gc()
  message("\nRunning Seurat")
  res.seurat <- run_Seurat(data, trueLabels = trueLabels)
  
  message("\n Generating Results")
  # curate the results
  
  log.rand <- adjustedRandIndex(trueLabels, obj@Labels$Seurat)
  multinom.rand <- adjustedRandIndex(trueLabels, obj@Labels$Scrna)
  negbinom.rand <- adjustedRandIndex(trueLabels, obj@Labels$SCTransform)
  multiview.rand <- adjustedRandIndex(trueLabels, obj@Labels$All)
  # Randomly.rand <- adjustedRandIndex(trueLabels, obj@Labels$Randomly)
  
  ensemble.mat <- cbind(obj@Labels$Seurat, obj@Labels$Scrna, obj@Labels$SCTransform)
  # ensemble.labs <- diceR::LCA(ensemble.mat)
  ensemble <- apply(ensemble.mat, MARGIN = 2, as.cl_partition)
  cons <- cl_consensus(as.cl_ensemble(ensemble), method = "SE", control = list(nruns = 50))
  ensemble.labs <- c(cl_class_ids(cons))
  obj@Labels$Ensemble <- ensemble.labs
  ensemle.rand <- adjustedRandIndex(trueLabels, ensemble.labs)
  
  res <- list(Log = log.rand, Multinom = multinom.rand, NegBinom = negbinom.rand, Ensemble = ensemle.rand, Multiview = multiview.rand,
                    Cidr = res.cidr$ARI, Monocle = res.monnocle$ARI, 
                    sscClust = res.sscClust$ARI, TSCAN = res.tscan$ARI, SC3 = res.sc3$ARI, Seurat = res.seurat$ARI)
  res <- as.data.frame(t(as.data.frame(res)))
  res <- setDT(res, keep.rownames = T)[]
  colnames(res) <- c("Method", "ARI")
  
  res$NumClusters <- c(Log = length(unique(obj@Labels$Seurat)), Multinom = length(unique(obj@Labels$Scrna)), NegBinom = length(unique(obj@Labels$SCTransform)),
                             Ensemble = length(unique(ensemble.labs)), Multiview = length(unique(obj@Labels$All)),
                             Cidr = length(unique(res.cidr$clusters)), Monocle = length(unique(res.monnocle$clusters)),
                             sscClust = length(unique(res.sscClust$clusters)), TSCAN = length(unique(res.tscan$clusters)), SC3 = length(unique(res.sc3$clusters)),
                             Seurat = length(unique(res.seurat$clusters)))
  
  res$NMI <- c(Log = NMI(trueLabels, obj@Labels$Seurat), Multinom = NMI(trueLabels, obj@Labels$Scrna), NegBinom = NMI(trueLabels, obj@Labels$SCTransform),
                     Ensemble = NMI(trueLabels, ensemble.labs), Multiview = NMI(trueLabels, obj@Labels$All),
                     Cidr = NMI(trueLabels, res.cidr$clusters), Monocle = NMI(trueLabels, res.monnocle$clusters),
                     sscClust = NMI(trueLabels, res.sscClust$clusters), TSCAN = NMI(trueLabels, res.tscan$clusters), SC3 = NMI(trueLabels, res.sc3$clusters),
                     Seurat = NMI(trueLabels, res.seurat$clusters))
  
  return(list(table = res, obj = obj))
}

# write a function to generate the plots
GeneratePlots <- function(obj, res.table, name, pca, k, res.dip){
  
  # plot heat map
  # res.table$Log[[1]] <- 0.001
  obj.heatmap <- pheatmap::pheatmap(as.matrix(res.dip), cluster_rows = F, cluster_cols = F, main = paste(name, "Dip heatmap"),
                                    display_numbers = T, number_color = "black")
  p.dip <- ggplotify::as.ggplot(obj.heatmap)

  # ARI
  p.ari <- res.table %>% 
    dplyr::mutate(Custom =  ifelse(Method == 'Log' | Method == 'NegBinom' | Method == 'Multinom' | Method == "Ensemble" | Method == "Multiview", T, F)) %>% 
    ggplot(aes(Method, ARI)) + geom_bar(stat="identity", aes(fill=Custom)) + labs(title = paste(name, "ARI")) + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                     axis.text.x = element_text(angle = 90))
  
  # NMI
  p.nmi <- res.table %>% 
    dplyr::mutate(Custom =  ifelse(Method == 'Log' | Method == 'NegBinom' | Method == 'Multinom' | Method == "Ensemble" | Method == "Multiview", T, F)) %>% 
    ggplot(aes(Method, NMI)) + geom_bar(stat="identity", aes(fill=Custom)) + labs(title = paste(name, "NMI")) + theme(plot.title = element_text(hjust = 0.5), 
                                                                                                                      axis.text.x = element_text(angle = 90))
  
  # Number of clusters
  g2 <- ggplot(res.table, aes(x=NumClusters, y=ARI, color=Method)) +
    theme_bw(base_size = 11, base_family = "") + guides(color=FALSE) +
    expand_limits(x=c(-5, max(res.table$NumClusters) + 5), y=c(0,1)) +
    scale_x_log10(breaks= seq(1, max(res.table$NumClusters) + 1, 3)) +
    geom_point(aes(NumClusters, ARI, color=Method)) +
    geom_text_repel(aes(NumClusters, ARI, label = Method, color=Method)) +
    geom_vline(xintercept=k, linetype = 2)
  
  # visualize the enseble
  evecs <- as.data.frame(obj@SCTransform$Tsne.Data)
  colnames(evecs) <- c("tSNE_1", "tSNE_2")
  Clusters <- as.factor(obj@Labels$Ensemble)
  p.evecs <- ggplot(evecs) + geom_point(aes(x = tSNE_1, y = tSNE_2, color = Clusters)) + 
    labs(title = paste(name, "Ensemble Clusters")) + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90))
  # combine the plots
  p.final <- plot_grid(pca, p.dip, p.ari, p.nmi, g2, p.evecs, labels = c("A", "B", "C", "D", "E", "F"))
  
  return(p.final)
  
}


EstimateK <- function(obj, method = "Seurat"){
  if(method == "Seurat"){
    d <- parallelDist(as.matrix(obj@Seurat$Princurve$lambda))
    modes <- modes::amps(d)
    k <- dim(modes$Peaks)[1]
  }
  
  if(method == "Scrna"){
    d <- parallelDist(as.matrix(obj@Scrna$Princurve$lambda))
    modes <- modes <- modes::amps(d)
    k <- dim(modes$Peaks)[1]
  }
  
  if(method == "SCTransform"){
    d <- parallelDist(as.matrix(obj@SCTransform$Princurve$lambda))
    modes <-  modes <- modes::amps(d)
    k <- dim(modes$Peaks)[1]
  }
  
  return(k)
}

PlotPCurve <- function(fit, method = "Log Model"){
  ord <- fit$tag
  df <- data.frame("Dim1" = fit$s[ord, 1], "Dim2" = fit$s[ord, 2])
  p.pcurve <- ggplot(data=df, aes(x=Dim1, y=Dim2)) +  geom_line(aes(color="red"), show.legend = FALSE) + geom_point(aes(fit$x[,1], fit$x[, 2], alpha = 1/10, color = trueLabels), show.legend = FALSE) + 
    ggtitle(paste(method, "Principal Curve")) + theme(plot.title = element_text(hjust = 0.5, size=15), axis.title=element_text(size=14))
  
  return(p.pcurve)
  
}

PlotPcDensity <- function(obj, test = "Dip"){
  Log <- as.data.frame(obj@Seurat$Sig.PCs.Data[, 1])
  Multinom <- as.data.frame(obj@Scrna$Sig.PCs.Data[, 1])
  NegBinom <- as.data.frame(obj@SCTransform$Sig.PCs.Data[, 1])
  
  df <- cbind(Log, Multinom, NegBinom)
  colnames(df) <- c("Log", "Multinom", "NegBinom")
  df.melted <- data.table::melt(df)
  colnames(df.melted) <- c("Model", "PC1")
  
  pvals <- as.data.frame(unlist(obj@Tests))
  if(test == "Dip"){
    pvals <-  pvals[rownames(pvals) %in% c("Seurat.Dip.PCA", "Scrna.Dip.PCA", "SCTransform.Dip.PCA"), ]
    p <- ggplot(df.melted, aes(PC1, colour = Model, fill = Model)) + geom_density(alpha = 0.2) +
      ggtitle("First Principal Component Density") + 
      theme(plot.title = element_text(hjust = 0.5, size=18), axis.title=element_text(size=18), legend.key.size = unit(1, "cm"), legend.title=element_text(size=12), legend.text=element_text(size=14))  +
      scale_color_discrete(name = "Dip Test Test P-value", labels = c(paste("Log:", round(pvals[1],2)), paste("Multinom:", round(pvals[2],2)),
                                                                       paste("NegBinom:", round(pvals[3],2))))
  } else if(test == "Silver"){
    pvals <- pvals[rownames(pvals) %in% c("Seurat.Silverman.PCA", "Scrna.Silverman.PCA", "SCTransform.Silverman.PCA"), ]
    p <- ggplot(df.melted, aes(PC1, colour = Model, fill = Model)) + geom_density(alpha = 0.2) +
      ggtitle("First Principal Component Density") + 
      theme(plot.title = element_text(hjust = 0.5, size=18), axis.title=element_text(size=18), legend.key.size = unit(1, "cm"), legend.title=element_text(size=12), legend.text=element_text(size=14))  +
      scale_color_discrete(name = "Silverman Test P-value", labels = c(paste("Log:", round(pvals[1],2)), paste("Multinom:", round(pvals[2],2)),
                                                                 paste("NegBinom:", round(pvals[3],2))))
  }
  
  return(p)
}


PlotPcurveDensity <- function(obj,test = "Dip"){
  Log <- as.data.frame(obj@Seurat$Princurve$lambda)
  Multinom <- as.data.frame(obj@Scrna$Princurve$lambda)
  NegBinom <- as.data.frame(obj@SCTransform$Princurve$lambda)
  
  df <- cbind(Log, Multinom, NegBinom)
  colnames(df) <- c("Log", "Multinom", "NegBinom")
  df.melted <- data.table::melt(df)
  colnames(df.melted) <- c("Model", "Lambda")
  
  # get the pvals
  pvals <- as.data.frame(unlist(obj@Tests))
  if(test == "Dip"){
   pvals <-  pvals[rownames(pvals) %in% c("Seurat.Dip.Princurve", "Scrna.Dip.Princurve", "SCTransform.Dip.Princurve"), ]
   # Generate the plot
   # Generate the plot
   p <- ggplot(df.melted, aes(Lambda, colour = Model, fill = Model)) + geom_density(alpha = 0.2) +
     ggtitle("Principal Curve Density") + 
     theme(plot.title = element_text(hjust = 0.5, size=18), axis.title=element_text(size=18), legend.key.size = unit(1, "cm"), legend.title=element_text(size=12), legend.text=element_text(size=14)) +
     scale_color_discrete(name = "Dip Test Test P-value", labels = c(paste("Log:", round(pvals[1],2)), paste("Multinom:", round(pvals[2],2)),
                                                                      paste("NegBinom:", round(pvals[3],2))))
  } else if(test == "Silver"){
    pvals <- pvals[rownames(pvals) %in% c("Seurat.Silverman.Princurve", "Scrna.Silverman.Princurve", "SCTransform.Silverman.Princurve"), ]
    # Generate the plot
    p <- ggplot(df.melted, aes(Lambda, colour = Model, fill = Model)) + geom_density(alpha = 0.2) +
      ggtitle("Principal Curve Density") + 
      theme(plot.title = element_text(hjust = 0.5, size=18), axis.title=element_text(size=18), legend.key.size = unit(1, "cm"), legend.title=element_text(size=12), legend.text=element_text(size=14)) +
      scale_color_discrete(name = "Silverman Test P-value", labels = c(paste("Log:", round(pvals[1],2)), paste("Multinom:", round(pvals[2],2)),
                                                                      paste("NegBinom:", round(pvals[3],2))))
  }
  

  
  return(p)
}

permutationPAMod <- function(dat,
                          B = 100, threshold = 0.05,
                          verbose = TRUE, seed = 1994, k = 10) {
  library(corpcor)
  if (!is.null(seed))
    set.seed(seed)
  n <- ncol(dat)
  m <- nrow(dat)
  
  uu <- RSpectra::svds(dat, k = k, nu = 0)
  ndf <- n - 1
  dstat <- uu$d[1:ndf]^2/sum(uu$d[1:ndf]^2)
  dstat0 <- matrix(0, nrow = B,
                   ncol = ndf)
  if (verbose == TRUE)
    message("Estimating a number of significant principal component: ")
  for (i in 1:B) {
    if (verbose == TRUE)
      cat(paste(i, " "))
    dat0 <- t(apply(dat, 1,
                    sample, replace = FALSE))
    uu0 <- RSpectra::svds(t(dat0), k = ndf, nu = 0)
    dstat0[i, ] <- uu0$d[1:ndf]^2/sum(uu0$d[1:ndf]^2)
  }
  p <- rep(1, n)
  for (i in 1:ndf) {
    p[i] <- mean(dstat0[, i] >=
                   dstat[i])
  }
  for (i in 2:ndf) {
    p[i] <- max(p[(i - 1)],
                p[i])
  }
  r <- sum(p < threshold)
  message(r)
  sigVecs <- irlba::prcomp_irlba(t(dat), n = max(3,r))$rotation
  
  
  return(list(r = r, p = p, sigVecs = sigVecs))
}



