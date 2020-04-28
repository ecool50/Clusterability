# seurat.genes <- rownames(obj@Seurat$Scaled.Data)
# scrna.genes <- rownames(obj@Scrna$Deviances)
# sct.genes <- rownames(obj@SCTransform$Scaled.Data)
# common.genes <- union(union(seurat.genes, scrna.genes), sct.genes)
# print(length(common.genes))


# 
# total.sigma <- seurat.sigma+scrna.sigma+sct.sigma

# cat("\nNow Generating Similarity Matrices")
# # generate the similarity matrices
# if(type == "All"){
#   cat("\nSEURAT")
#   # d1 <- obj@Seurat$Umap.Data
#   # Seurat
#   # d1 <- parallelDist(t(obj@Seurat$Scaled.Data))
#   # d1 <- abs(coop::pcor((obj@Seurat$Scaled.Data)))
#   d1 <- coop::covar((obj@Seurat$Scaled.Data))
#   
#   # d1 <- .rbfkernel_b((obj@Seurat$Scaled.Data), K = 15)
#   # d1 <- .rbfkernel_b(as.matrix(d1))
#   # d1 <- parallelDist(obj@Seurat$Sig.PCs.Data)
#   # d2 <- parallelDist(obj@Seurat$Umap.Data)
#   # d4 <- 1 - WGCNA::cor1(t(obj@Seurat$Sig.PCs.Data))
#   # d7 <- 1 - WGCNA::cor1((obj@Seurat$Scaled.Data))
#   # d1 <- SC3::norm_laplacian(as.matrix(d1))
# 
#   # Scrna
#   cat("\nSCRNA")
#   # d2 <- parallelDist(t(obj@Scrna$Deviances))
#   # d2 <- abs(coop::pcor((obj@Scrna$Deviances)))
#   d2 <- coop::covar((obj@Scrna$Deviances))
#   rownames(obj@Scrna$Umap.Data) <- colnames(obj@Scrna$Deviances)
#   rownames(obj@Scrna$Sig.PCs.Data) <- colnames(obj@Scrna$Deviances)
#   # d2 <-  .rbfkernel_b((obj@Scrna$Deviances), K = 15)
#   # d2<- .LaplacianNg(as.matrix(d2))
#     # .rbfkernel_b(obj@Scrna$Deviances, scrna.sigma)
#   # d4 <- parallelDist(obj@Scrna$Sig.PCs.Data)
#   # d5 <- parallelDist(obj@Scrna$Umap.Data)
#   # d2 <- obj@Scrna$Umap.Data
#   # d5 <- 1 - WGCNA::cor1(t(obj@Scrna$Sig.PCs.Data))
#   # d8 <- 1 - WGCNA::cor1((obj@Scrna$Scaled.Data))
#   # d2 <- SC3::norm_laplacian(as.matrix(d2))
# 
#   # SCTransform
#   cat("\nSCTRANSFORM")
#   # d3 <- parallelDist(t(obj@SCTransform$Scaled.Data))
#   # d3 <- abs(coop::pcor((obj@SCTransform$Scaled.Data)))
#   # d3 <- .rbfkernel_b((obj@SCTransform$Scaled.Data), K = 15)
#   # d3 <- .LaplacianNg(as.matrix(d3))
#   d3 <- coop::covar((obj@SCTransform$Scaled.Data))
#   # d7 <- parallelDist(obj@SCTransform$Sig.PCs.Data)
#   # d8 <- parallelDist(obj@SCTransform$Umap.Data)
#   # d3 <- obj@SCTransform$Umap.Data
#   # d6 <- 1 - WGCNA::cor1(t(obj@SCTransform$Sig.PCs.Data))
#   # d8 <- 1 - WGCNA::cor1((obj@SCTransform$Scaled.Data))
#   # d3 <- SC3::norm_laplacian(as.matrix(d3))
#   
#   ng <- c(nrow(d1), nrow(d2), nrow(d3))
#   mat <- list(d1,d2,d3)
#   
#   res <- .prepMatrix(covs = mat, ng = ng, ncomp = ncomp)
#   cat("\nNow Running MultiView")
#   # res <- .integrate_similarity_matrices(mat, KNNs_p = 15, diffusion_iters = 4)
#   # res <- .LaplacianNg(res)
#   x <- as.matrix(as.data.frame(res$v0))
#   
#   # x <- .cpca_stepwise_eigen(mat, ng, ncomp = ncomp, start = "eigen", symmetric = TRUE, maxit = 10, tol = 1e-4, verbose = 2)$CPC
#   
#   
#   # x <- eigenPower(A = res$A, v0 = res$v0, verbose = 2, sparse = T, ncomp = ncomp)$vectors
#   # x <- as.matrix(.SigPCs(res$A, seed = 1994, ndf = 50)$sigVecs)
#   # x <- multiview::mvsc(mat, ncomp)$evectors
#   
# 
# } else if(type == "Seurat"){
#     cat("\nNow Running Eigen Decomposition")
#     # d1 <- .ng_kernel((obj@Seurat$Scaled.Data))
#       # d1 <- coop::covar((obj@Seurat$Scaled.Data))
#     x <- obj@Seurat$Sig.PCs.Data[, 1:5]
#     # x <- .estkTW(obj@Seurat$Scaled.Data, ncomp = ncomp)
#       # x <- RSpectra::eigs_sym(d1, k = ncomp)$vectors
#   # x <- as.matrix(.SigPCs(d1, seed = 1994, ndf = 25)$sigVecs)
#   # x <- t(obj@Seurat$Scaled.Data)
# } else if(type == "Scrna"){
#     cat("\nNow Running Eigen Decomposition")
#   rownames(obj@Scrna$Sig.PCs.Data) <- colnames(obj@Scrna$Deviances)
#   # d2 <- .ng_kernel((obj@Scrna$Deviances))
#  # d2 <- coop::covar(scale(obj@Scrna$Deviances))
#   x <- obj@Scrna$Sig.PCs.Data[, 1:5]
#   # x <- .estkTW(obj@Scrna$Deviances, ncomp = ncomp)
#  # x <- RSpectra::eigs_sym(d2, k = ncomp)$vectors
#   # x <- as.matrix(.SigPCs(t(obj@Scrna$Deviances), seed = 1994, ndf = 25)$sigVecs)
#       # RSpectra::eigs_sym(as.matrix(cov_mat), k = ncomp)$vectors
# } else if(type == "SCTransform"){
#   cat("\nNow Running Eigen Decomposition")
#   # d3 <- .ng_kernel((obj@SCTransform$Scaled.Data))
#  # d3 <- coop::covar((obj@SCTransform$Scaled.Data))
#  # x <- RSpectra::eigs_sym(d3, k = ncomp)$vectors
#   x <- obj@SCTransform$Sig.PCs.Data[, 1:5]
#   # x <- .estkTW(obj@SCTransform$Scaled.Data, ncomp = ncomp)
#  # x <- as.matrix(.SigPCs(t(obj@SCTransform$Scaled.Data), seed = 1994, ndf = 25)$sigVecs)
#   
# }

# compute the common eigenvectors
# mat <- array(dim = c(n,n,3))
# mat[,,1] <- as.matrix(d1)
# mat[,,2] <- as.matrix(d2)
# mat[,,3] <- as.matrix(d3)
# mat[,,4] <- as.matrix(d4)
# mat[,,5] <- as.matrix(d5)
# mat[,,6] <- as.matrix(d6)
# mat[,,7] <- as.matrix(d7)
# mat[,,8] <- as.matrix(d8)
# # mat[,,9] <- as.matrix(d9)
# 






# if(type == "All"){
#   obj@Tests$Multiview$All <-  list(in_mat = x,
#                                     nd_type = nd_type,
#                                     p_emp = p_emp,
#                                     idx_hc = idx_hc,
#                                     hc_dat = hc_dat,
#                                     pd_map = pd_map,
#                                     evecs = x)
# } else if(type == "Seurat"){
#    obj@Tests$Multiview$Seurat <-  list(in_mat = x,
#                                     nd_type = nd_type,
#                                     p_emp = p_emp,
#                                     idx_hc = idx_hc,
#                                     hc_dat = hc_dat,
#                                     pd_map = pd_map,
#                                     evecs = x)
# } else if(type == "Scrna"){
#   obj@Tests$Multiview$Scrna <-  list(in_mat = x,
#                                    nd_type = nd_type,
#                                    p_emp = p_emp,
#                                    idx_hc = idx_hc,
#                                    hc_dat = hc_dat,
#                                    pd_map = pd_map,
#                                    evecs = x)
# } else if(type == "SCTransform"){
#   obj@Tests$Multiview$SCTransform <-  list(in_mat = x,
#                                      nd_type = nd_type,
#                                      p_emp = p_emp,
#                                      idx_hc = idx_hc,
#                                      hc_dat = hc_dat,
#                                      pd_map = pd_map,
#                                      evecs = x)
# }
