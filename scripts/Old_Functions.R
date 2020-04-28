# function to run a significant test on the result of consensus clustering
.Consensus <- function(obj, method = "Seurat"){
  if(!(method %in% c("Scrna", "Seurat", "SCTransform"))){
    stop(paste(type, "is not a valid option"))
  }
  
  if(method == "Seurat"){
    # compute consensus for the scaled data
    con_output <- .Boot_Con(t(obj@Seurat$Scaled.Data), obj, distance="euclidean", method="ward",nboot=100)
    # con_dend <- as.dendrogram(con_output$dendrogram)
    # get the consensus labels
    assignments <- cutree(con_output$dendrogram, 2)
    # run the statistical significance test
    res_raw <- sigclust2::sigclust(t(obj@Seurat$Scaled.Data), labels = assignments)$p_norm
    
    # compute the consensus for the sig pc data
    con_output <- .Boot_Con(obj@Seurat$Sig.PCs.Data, obj, distance = "euclidean", method = "ward", nboot = 100)
    # con_dend <- as.dendrogram(con_output$dendrogram)
    # get the consensus labels
    assignments <- cutree(con_output$dendrogram, 2)
    # run the statistical significance test
    res_pc <- sigclust2::sigclust(obj@Seurat$Sig.PCs.Data, labels = assignments)$p_norm
    # update the modality object
    obj@Tests$Seurat$Consensus <- list("Raw" = round(res_raw,digits = 3), "PC" = round(res_pc,digits = 3))
    
  }
  
  if(method == "Scrna"){
    # compute consensus for the scaled data
    con_output <- .Boot_Con(t(obj@Scrna$Deviances), obj, distance="euclidean", method="ward",nboot=100)
    # con_dend <- as.dendrogram(con_output$dendrogram)
    # get the consensus labels
    assignments <- cutree(con_output$dendrogram, 2)
    # run the statistical significance test
    res_raw <- sigclust2::sigclust(t(obj@Scrna$Deviances), labels = assignments)$p_norm
    
    # compute the consensus for the sig pc data
    con_output <- .Boot_Con(obj@Scrna$Sig.PCs.Data, obj, distance = "euclidean", method = "ward", nboot = 100)
    # con_dend <- as.dendrogram(con_output$dendrogram)
    # get the consensus labels
    assignments <- cutree(con_output$dendrogram, 2)
    # run the statistical significance test
    res_pc <- sigclust2::sigclust(obj@Seurat$Sig.PCs.Data, labels = assignments)$p_norm
    # update the modality object
    obj@Tests$Scrna$Consensus <- list("Raw" = round(res_raw,digits = 3), "PC" = round(res_pc,digits = 3))
    
  }
  
  
  if(method == "SCTransform"){
    # compute consensus for the scaled data
    con_output <- .Boot_Con(t(obj@SCTransform$Scaled.Data), obj, distance="euclidean", method="ward",nboot=100)
    # con_dend <- as.dendrogram(con_output$dendrogram)
    # get the consensus labels
    assignments <- cutree(con_output$dendrogram, 2)
    # run the statistical significance test
    res_raw <- sigclust2::sigclust(t(obj@SCTransform$Scaled.Data), labels = assignments)$p_norm
    
    # compute the consensus for the sig pc data
    con_output <- .Boot_Con(obj@SCTransform$Sig.PCs.Data, obj, distance = "euclidean", method = "ward", nboot = 100)
    # con_dend <- as.dendrogram(con_output$dendrogram)
    # get the consensus labels
    assignments <- cutree(con_output$dendrogram, 2)
    # run the statistical significance test
    res_pc <- sigclust2::sigclust(obj@SCTransform$Sig.PCs.Data, labels = assignments)$p_norm
    # update the modality object
    obj@Tests$SCTransform$Consensus <- list("Raw" = round(res_raw,digits = 3), "PC" = round(res_pc,digits = 3))
  }
  
  # Return the object
  return(obj)
}


# Run clustering and generate randIndex values
.ClusterDat <- function(obj, sce, trueLabels, sc3_num_clust, All = F){
  
  # Seurat
  cat("\nNow Running Seurat")
  seurat.obj <- obj@Seurat$obj
  seurat.obj <- Seurat::FindNeighbors(seurat.obj, dims = 1:ncomp)
  seurat.obj <- Seurat::FindClusters(seurat.obj, resolution = 0.8, random.seed = 1994)
  res.seurat <- Seurat::Idents(seurat.obj) 
  
  
  # run clustering
  cat("\nNow Running MultiView")
  
  # everything else
  res.sc3 <- as.numeric(colData(sce)[, paste0("sc3_", sc3_num_clust, "_clusters")])
  cat("\nNow Running CIDR")
  res.cidr <- run_CIDR(as.matrix(obj@Data$cleaned), ncomp = ncomp)
  cat("\nNow Running RaceID3")
  res.race <- run_RaceID3(as.matrix(obj@Data$cleaned))
  cat("\nNow Running TSCAN")
  res.tscan <- run_TSCAN(as.matrix(obj@Data$cleaned))
  cat("\nNow Running Ascend")
  res.ascend <- run_ascend(counts(sce))
  cat("\nNow Running Monocle")
  res.monocle <- run_monocle3(as.matrix(obj@Data$cleaned))
  cat("\nNow Running Sincell")
  res.sincell <- run_sincell(as.matrix(obj@Data$cleaned))
  # cat("\nNow Running sscClust")
  # res.sscClust <- run_sscClust(as.matrix(obj@Data$cleaned))
  # res.flosom <- FlowSOM::metaClustering_consensus(fSOM$map$codes, k = numNonNoiseClusters)
  # res.flosom <- res.flosom[fSOM$map$mapping[, 1]]
  # res.rtsneKmeans <- kmeans(obj@Scrna$Tsne.Data, centers = numNonNoiseClusters, nstart = 25)$cluster
  # res.umapKmeans <- kmeans(obj@Scrna$Umap.Data, centers = numNonNoiseClusters, nstart = 25)$cluster
  #res.pcaReduce <- c(clue::cl_class_ids(cons))
  # res.simlr <- res.simlr$y$cluster
  # res.monocle <- clusterCells(cds, num_clusters = numNonNoiseClusters + 1, method = "densityPeak")$Cluster
  
  gc()
  res.sc3.rand <- mclust::adjustedRandIndex(trueLabels, res.sc3)
  res.cidr.rand <- mclust::adjustedRandIndex(trueLabels, res.cidr)
  res.seurat.rand <- mclust::adjustedRandIndex(trueLabels, res.seurat)
  res.ascend.rand <- mclust::adjustedRandIndex(trueLabels, unlist(res.ascend))
  res.tscan.rand <- mclust::adjustedRandIndex(trueLabels, res.tscan)
  # res.simlr.rand <- adjustedRandIndex(trueLabels, res.simlr)
  res.race.rand <- adjustedRandIndex(trueLabels, unlist(res.race))
  res.monocle.rand <- adjustedRandIndex(trueLabels, res.monocle)
  res.sincell.rand <- adjustedRandIndex(trueLabels, res.sincell)
  # res.sscClust.rand <- adjustedRandIndex(trueLabels, res.sscClust)
  
  
  if(All == T){
    Seurat.Sig <- adjustedRandIndex(obj@Labels$Seurat, trueLabels)
    Scrna.Sig <- adjustedRandIndex(obj@Labels$Scrna, trueLabels)
    SCTransform.Sig <- adjustedRandIndex(obj@Labels$SCTransform, trueLabels)
    All.Sig <- adjustedRandIndex(obj@Labels$All, trueLabels)
    Ensemble.Sig <- adjustedRandIndex(diceR::LCA(cbind(obj@Labels$Seurat, obj@Labels$Scrna, obj@Labels$SCTransform, obj@Labels$All)),
                                      trueLabels)
    
    rand.res <- list(Sc3 = res.sc3.rand, Multiview.Sig = All.Sig, Seurat = res.seurat.rand, Seurat.Sig = Seurat.Sig, 
                     Scrna.Sig = Scrna.Sig, SCTransform.Sig = SCTransform.Sig, Ensemble.Sig = Ensemble.Sig, RaceID3 = res.race.rand,
                     TSCAN = res.tscan.rand, Ascend = res.ascend.rand, CIDR = res.cidr.rand, Monocle = res.monocle.rand,
                     Sincell = res.sincell.rand
    )
  } else{
    rand.res <- list(Sc3 = res.sc3.rand, Multiview.Sig = res.multi.sig, Seurat = res.seurat.rand, CIDR = res.cidr.rand,
                     RaceID3 = res.race.rand, TSCAN = res.tscan.rand, Ascend = res.ascend.rand, Monocle = res.monocle.rand,
                     Sincell = res.sincell.rand)
  }
  
  return(rand.res)
}

ComputeProportions <- function(obj){
  
  # compute a result table for each of the test types
  # Seurat
  # seurat_raw_dips <- as.data.frame(unlist(obj@Tests$Seurat$Raw))
  # colnames(seurat_raw_dips) <- c("Dip.Pvalue")
  seurat_pc_dips <- as.data.frame(unlist(obj@Tests$Seurat$PC))
  colnames(seurat_pc_dips) <- c("Dip.Pvalue")
  seurat_mdh_dip <- as.data.frame(unlist(obj@Tests$Seurat$MDH))
  colnames(seurat_mdh_dips) <- c("Dip.Pvalue")
  
  # Scrna
  scrna_raw_dips <- as.data.frame(unlist(obj@Tests$Scrna$Raw))
  colnames(scrna_raw_dips) <- c("Dip.Pvalue")
  scrna_pc_dips <- as.data.frame(unlist(obj@Tests$Scrna$PC))
  colnames(scrna_pc_dips) <- c("Dip.Pvalue")
  scrna_mdh_dips <- as.data.frame(unlist(obj@Tests$Scrna$MDH))
  colnames(scrna_mdh_dips) <- c("Dip.Pvalue")
  
  # compute proportions
  # Seurat
  seurat_raw <- count(seurat_raw_dips$Dip.Pvalue < 0.05)/length(seurat_raw_dips$Dip.Pvalue)
  seurat_pc <-   count(seurat_pc_dips$Dip.Pvalue < 0.05)/length(seurat_pc_dips$Dip.Pvalue)
  seurat_mdh <-  count(seurat_mdh_dips$Dip.Pvalue < 0.05)/length(seurat_mdh_dips$Dip.Pvalue)
  seurat_prop <- data.frame(seurat_raw, seurat_pc, seurat_mdh, row.names = c("Proportion"))
  obj@Tests$Seurat$Proportions <- seurat_prop
  
  
  # Scrna
  scrna_raw <- count(scrna_raw_dips$Dip.Pvalue < 0.05)/length(scrna_raw_dips$Dip.Pvalue)
  scrna_pc <-   count(scrna_pc_dips$Dip.Pvalue < 0.05)/length(scrna_pc_dips$Dip.Pvalue)
  scrna_mdh <-  count(scrna_mdh_dips$Dip.Pvalue < 0.05)/length(scrna_mdh_dips$Dip.Pvalue)
  scrna_prop <- data.frame(scrna_raw, scrna_pc, scrna_mdh, row.names = c("Proportion"))
  obj@Tests$Scrna$Proportions <- scrna_prop
  
  # return the object
  return(obj)
  
}

#' Estimate the optimal k for k-means clustering
#' 
#' The function finds the eigenvalues of the sample covariance matrix. 
#' It will then return the number of significant eigenvalues according to 
#' the Tracy-Widom test.
#' 
#' @param dataset processed input expression matrix.
#' @return an estimated number of clusters k
estkTW <- function(dataset) {
  
  p <- ncol(dataset)
  n <- nrow(dataset)
  
  # compute Tracy-Widom bound
  #x <- scale(dataset)
  muTW <- (sqrt(n - 1) + sqrt(p))^2
  sigmaTW <- (sqrt(n - 1) + sqrt(p)) * (1/sqrt(n - 1) + 1/sqrt(p))^(1/3)
  sigmaHatNaive <- Rfast::mat.mult(Rfast::transpose(dataset), dataset)  # x left-multiplied by its transpose
  bd <- 3.273 * sigmaTW + muTW  # 3.2730 is the p=0.001 percentile point for the Tracy-Widom distribution
  
  # compute eigenvalues and return the amount which falls above the bound
  evals <- eigen(sigmaHatNaive, symmetric = TRUE, only.values = TRUE)$values
  k <- 0
  for (i in 1:length(evals)) {
    if (evals[i] > bd) {
      k <- k + 1
    }
  }
  return(k)
}

# write a function to shuffle a matrix. This function has been borrowed from the Seurat package.

.MatrixRowShuffle <- function(x) {
  x2 <- x
  x2 <- t(x = x)
  ind <- order(c(col(x = x2)), runif(n = length(x = x2)))
  x2 <- matrix(
    data = x2[ind],
    nrow = nrow(x = x),
    ncol = ncol(x = x),
    byrow = TRUE
  )
  return(x2)
}


Dim_Reduction <- function(data, nperm = 1000, nPCs = 25, seed = 1994, prop = 0.1){
  # pca_out<- prcomp(data, scale. = F)
  eigenperm <- data.frame(matrix(NA, nperm, nPCs))
  eigenreal <- data.frame(matrix(NA, nperm, nPCs))
  n<- ncol(data)
  #data_i<- data.frame(matrix(NA, nrow(data), ncol(data)))
  #data_mod <- t(data)
  for (j in 1: nperm){
    random_genes <- sample(x = colnames(x = data),
                           size = ceiling(ncol(x = data) * prop))
    data_mod <- t(data)[random_genes, ]
    
    data_mod.shuff <- .MatrixRowShuffle(data_mod) 
    
    pca.perm <- prcomp_irlba(t(data_mod.shuff), n = nPCs, retx = TRUE, center = F, scale. = FALSE, fastpath = TRUE)
    pca.real <- prcomp_irlba(t(data_mod), n = nPCs, retx = TRUE, center = F, scale. = FALSE, fastpath = TRUE)
    # prcomp(t(data_mod), scale. = F, rank. = nPCs)
    eigenperm[j,] <- pca.perm$sdev^2
    eigenreal[j, ] <- pca.real$sdev^2
  }
  
  # compute the principal components on the full dataset
  pvals <- colSums(eigenperm > eigenreal)/length(eigenperm)
  
  data_pc <- prcomp_irlba(data, n = nPCs, retx = TRUE, center = F, scale. = FALSE, fastpath = TRUE)
  #   # prcomp(data, scale. = F,  rank. = nPCs)
  # res <- data.frame(Random_Eigenvalues = sapply(eigenperm, quantile, 0.95))
  # res$PCs <- colnames(data_pc$rotation)
  # res$Eigenvalues <- data_pc$sdev[1:nPCs]^2
  
  # get indices of significant pcs
  sig_indices <- which(pvals < 0.05)
  num_sig_pcs <- length(sig_indices)
  
  cat(paste("\nNumber of significant PCs:", num_sig_pcs))
  # if we have more significant PCS than required
  if(nPCs >= num_sig_pcs && num_sig_pcs >= 10){
    # run PCA and only return the number specified
    cat("\nNumber of significant PCs is at least the number provided. Returning the significant PCs")
    pca_res <- data_pc$x[, sig_indices]
  }
  
  # if we have less PCs than required
  if(nPCs < num_sig_pcs){
    cat("\nNumber of significant PCs is greater than the number provided. Returning the number of PCs specified.")
    pca_res <- data_pc$x[, sig_indices[1:nPCs]]
  }
  
  if(num_sig_pcs < 10){
    cat("\nNumber of significant PCs is less than 10. Returning the top 10 PCs.")
    pca_res <- data_pc$x[, 1:15]
  }
  
  # run MDH on the result
  sol <- mddr(data, 2)
  
  
  # return the pca object
  return(list("Sig.PCs" = pca_res, "MDH" = sol, "pvals" = pvals))
  
}


# create a function to the silverman test
.RunSilver <- function(dat, type = "PC"){
  if(!(type %in% c("PC", "Raw", "MDH", "Tsne", "Umap", "Multiview"))){
    stop(paste(type, "is not a valid option"))
  }
  
  # if we are dealing with PC data
  if(type == "PC"){
    cat("\nRunning tests on PC data\n")
    # compute the euclidean distances
    eucl_dist <- as.matrix(parallelDist::parallelDist(dat))
    eucl_dist <- t(eucl_dist)[lower.tri(t(eucl_dist))]
    eucl_silver <- multimode::modetest(eucl_dist, method = "HY", B = 10)$p.value
    rm(eucl_dist)
    gc()
    
    # compute the manhattan distances
    man_dist <- as.matrix(parallelDist::parallelDist(dat, method = "manhattan"))
    man_dist <- t(man_dist)[lower.tri(t(man_dist))]
    man_silver <- multimode::modetest(man_dist, method = "HY", B = 10)$p.value
    rm(man_dist)
    gc()
    
    # compute the pearson correlation distances
    pear_corr <- as.matrix(1 - WGCNA::cor1(t(dat), verbose = F))
    pear_corr <- t(pear_corr)[lower.tri(t(pear_corr))]
    pear_silver <- multimode::modetest(pear_corr, method = "HY", B = 10)$p.value
    rm(pear_corr)
    gc()
    
    # run dip test on the first component
    # comp.1 <- diptest::dip.test(dat[, 1])$p.value
    
    # return the results
    return(list("Euclidean" = eucl_silver, "Manhattan" = man_silver, "Pearson" = pear_silver))
  }
  
  # if we are computing on the raw data
  if(type == "Raw"){
    cat("\nRunning tests on raw data\n")
    # compute the euclidean distances
    eucl_dist <- as.matrix(parallelDist::parallelDist(t(dat)))
    eucl_dist <- t(eucl_dist)[lower.tri(t(eucl_dist))]
    eucl_silver <- multimode::modetest(eucl_dist, method = "HY", B = 10)$p.value
    rm(eucl_dist)
    gc()
    
    # compute the manhattan distances
    man_dist <- as.matrix(parallelDist::parallelDist(t(dat), method = "manhattan"))
    man_dist <- t(man_dist)[lower.tri(t(man_dist))]
    man_silver <- multimode::modetest(man_dist, method = "HY", B = 10)$p.value
    rm(man_dist)
    gc()
    
    # compute the pearson correlation distances
    pear_corr <- as.matrix(1 - WGCNA::cor1(dat, verbose = F))
    pear_corr <- t(pear_corr)[lower.tri(t(pear_corr))]
    pear_silver <- multimode::modetest(pear_corr, method = "HY", B = 10)$p.value
    rm(pear_corr)
    gc()
    
    return(list("Euclidean" = eucl_silver, "Manhattan" = man_silver, "Pearson" = pear_silver))
  }
  
  # if we are using the projected data
  if(type == "MDH"){
    cat("\nRunning tests on MDH data\n")
    # compute the euclidean distances
    eucl_dist <- as.matrix(parallelDist::parallelDist(dat))
    eucl_dist <- t(eucl_dist)[lower.tri(t(eucl_dist))]
    eucl_silver <- multimode::modetest(eucl_dist, method = "HY", B = 10)$p.value
    rm(eucl_dist)
    gc()
    
    # compute the manhattan distances
    man_dist <- as.matrix(parallelDist::parallelDist(dat, method = "manhattan"))
    man_dist <- t(man_dist)[lower.tri(t(man_dist))]
    man_silver <- multimode::modetest(man_dist, method = "HY", B = 10)$p.value
    rm(man_dist)
    gc()
    
    # compute the pearson correlation distances
    # pear_corr <- 1 - WGCNA::cor1(t(dat), verbose = F)
    # pear_dip <- diptest::dip.test(pear_corr)$p.value
    # rm(pear_corr)
    gc()
    # return the results
    return(list("Euclidean" = eucl_silver, "Manhattan" = man_silver))
    
  }
  
  # if we are using Tsne
  if(type == "Tsne"){
    cat("\nRunning tests on Tsne data\n")
    # compute the euclidean distances
    eucl_dist <- parallelDist::parallelDist(dat)
    eucl_silver <- silvermantest::silverman.test(eucl_dist, k = 1, R = 50, adjust = T)@p_value
    rm(eucl_dist)
    gc()
    
    # compute the manhattan distances
    man_dist <- parallelDist::parallelDist(dat, method = "manhattan")
    man_silver <- silvermantest::silverman.test(man_dist, k = 1, R = 50, adjust = T)@p_value
    rm(man_dist)
    gc()
    
    # run dip test on the first component
    comp.1 <- silvermantest::silverman.test(dat[, 1], k = 1, R = 50, adjust = T)@p_value
    # return the results
    return(list("Euclidean" = eucl_silver, "Manhattan" = man_silver, "Comp_1" = comp.1))
    
  }
  
  # if we are using UMAP
  if(type == "Umap"){
    cat("\nRunning tests on Umap data\n")
    # compute the euclidean distances
    eucl_dist <- as.matrix(parallelDist::parallelDist(dat))
    eucl_dist <- t(eucl_dist)[lower.tri(t(eucl_dist))]
    eucl_silver <- multimode::modetest(eucl_dist, method = "HY", B = 10)$p.value
    rm(eucl_dist)
    gc()
    
    # compute the manhattan distances
    man_dist <- as.matrix(parallelDist::parallelDist(dat, method = "manhattan"))
    man_dist <- t(man_dist)[lower.tri(t(man_dist))]
    man_silver <- multimode::modetest(man_dist, method = "HY", B = 10)$p.value
    rm(man_dist)
    gc()
    
    # run dip test on the first component
    # comp.1 <- diptest::dip.test(dat[, 1])$p.value
    
    # return the results
    return(list("Euclidean" = eucl_silver, "Manhattan" = man_silver))
  }
  
  if(type == "Multiview"){
    cat("\nRunning tests on Multiview data\n")
    # compute the euclidean distances
    eucl_dist <- parallelDist(dat)
    eucl_dip <- dip.test(eucl_dist)$p.value
    rm(eucl_dist)
    gc()
    
    # compute the manhattan distances
    man_dist <- parallelDist(dat, method = "manhattan")
    man_dip <- dip.test(man_dist)$p.value
    rm(man_dist)
    gc()
    
    # compute the pearson correlation distances
    # pear_corr <- 1 - WGCNA::cor1(t(dat), verbose = F)
    # pear_dip <- diptest::dip.test(pear_corr)$p.value
    # rm(pear_corr)
    # gc()
    # return the results
    return(list("Euclidean" = eucl_dip, "Manhattan" = man_dip))
  }
}



.RunHopkins <- function(dat, type = "PC"){
  if(!(type %in% c("PC", "Raw", "MDH", "Tsne", "Umap", "Multiview"))){
    stop(paste(type, "is not a valid option"))
  }
  
  if(type == "PC"){
    res <- factoextra::get_clust_tendency(dat, n = 100, graph = F)
  }
  
  if(type == "Raw"){
    res <- factoextra::get_clust_tendency(t(dat), n = 100, graph = F)
  }
  
  if(type == "MDH"){
    res <- factoextra::get_clust_tendency(dat, n = 100, graph = F)
  }
  
  if(type == "Tsne"){
    res <- factoextra::get_clust_tendency(dat, n = 100, graph = F)
  }
  
  if(type == "Umap"){
    res <- factoextra::get_clust_tendency(dat, n = 100, graph = F)
  }
  
  return(list("Hopkins" = res$hopkins_stat))
}



.RunMultiview <- function(obect, method = "Seurat", ncomp = 5){
  if(!(method %in% c("Seurat", "Scrna", "SCTransform"))){
    stop(paste(type, "is not a valid option"))
  }
  
  if(method == "Seurat"){
    cat("Running Multiview on Seurat")
    # extract the data
    x1 <- obj@Seurat$MDH$fitted
    x2 <- t(obj@Seurat$Scaled.Data)
    x3 <- obj@Seurat$Sig.PCs.Data
    x4 <- obj@Seurat$Tsne.Data
    x5 <- obj@Seurat$Umap.Data
    
    # run multiview
    res.multi <- mvsc(list(x1,x3,x4,x5), ncomp)$evectors
    
    # return the results
    obj@Seurat$Multiview.Data <- res.multi
  }
  
  if(method == "Scrna"){
    cat("Running Multiview on Scrna")
    # extract the data
    x1 <- obj@Scrna$MDH$fitted
    x2 <- t(obj@Scrna$Deviances)
    x3 <- obj@Scrna$Sig.PCs.Data
    x4 <- obj@Scrna$Tsne.Data
    x5 <- obj@Scrna$Umap.Data
    
    # run multiview
    res.multi <- mvsc(list(x1,x3,x4,x5), ncomp)$evectors
    
    # return the results
    obj@Scrna$Multiview.Data <- res.multi
  }
  
  if(method == "SCTransform"){
    cat("Running Multiview on SCTransform")
    # extract the data
    x1 <- obj@SCTransform$MDH$fitted
    x2 <- t(obj@SCTransform$Scaled.Data)
    x3 <- obj@SCTransform$Sig.PCs.Data
    x4 <- obj@SCTransform$Tsne.Data
    x5 <- obj@SCTransform$Umap.Data
    
    # run multiview
    res.multi <- mvsc(list(x1,x3,x4,x5), ncomp)$evectors
    
    # return the results
    obj@SCTransform$Multiview.Data <- res.multi
  }
  
  # return the results
  return(obj)
  
}


# Create a function to do dimensionality reduction and returns the significant PCS using the jackstraw method of permuation.

# Dim_Reduction <- function(dat, nPCs = 25, seed = 1994){
#   
#   # if(nrow(dat) <= ncol(dat)){
#   #   num_sig_pcs <- permutationPA(t(dat), B = 100, threshold = 0.05, verbose = TRUE,
#   #                                seed = seed)$r
#   # }
#   # else{
#   #   # compute the number of significant PCS using a 100 permutations
#   #   num_sig_pcs <- permutationPA(dat, B = 100, threshold = 0.05, verbose = TRUE,
#   #                                seed = seed)$r
#   # }
#   
#   # run sample pca
#   spca <- SamplePCA(dat, usecor = F, center = F)
#   
#   # Compute the significant PCs
#   spca_ag <- AuerGervini(spca)
#   num_sig_pcs <- agDimension(spca_ag, agfun = agDimSpectral)
#   cat(paste("Number of significant PCs:", num_sig_pcs))
#   # if we have more significant PCS than required
#   if(nPCs >= num_sig_pcs && num_sig_pcs >= 10){
#     # run PCA and only return the number specified
#     cat("\nNumber of significant PCs is atleast the number provided. Returning the significant PCs")
#     pca_res <- spca@scores[, 1:num_sig_pcs]
#   }
#   
#   # if we have less PCs than required
#   if(nPCs <= num_sig_pcs){
#     cat("\nNumber of significant PCs is greater than the number provided. Returning the number of PCs specified.")
#     pca_res <- spca@scores[, 1:nPCs]
#   }
#   
#   if(num_sig_pcs < 10){
#     cat("\nNumber of significant PCs is less than 10. Returning the top 10 PCs.")
#     pca_res <- spca@scores[, 1:15]
#   }
#   
#   # run MDH on the result
#   if(num_sig_pcs >= 10){
#     sol <- mddr(t(dat), num_sig_pcs, v0 = spca@components[, 1])
#   }
#   else{
#     sol <- mddr(t(dat), 10, v0 = spca@components[, 1])
#   }
#   
#   # return the pca object
#   return(list("Sig.PCs" = pca_res, "MDH" = sol))
# }
