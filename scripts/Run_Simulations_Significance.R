# Simulation Setup 1
# Balanced
Ensemble <- c()
Multiview <- c()
EnsembleNumClust <- c()
MultiviewNumClust <- c()
for(i in seq(5000, 15000, 5000)){
  for(j in c(4, 8, 16)){
    message(paste("Now running for numCells:", i))
    message(paste("Now running for numClusters:", j))
    sce <- Setup_1(numCells = i, numClusters = j)
    meta.data <- colData(sce)
    trueLabels <- as.numeric(as.factor(meta.data$Group))
    
    # Run singnificance testing
    res <- RunSignificanceTest(counts(sce))
    
    # compute results
    ensemble.mat <- cbind(res@Labels$Seurat, res@Labels$Scrna, res@Labels$SCTransform)
    ensemble <- apply(ensemble.mat, MARGIN = 2, as.cl_partition)
    cons <- cl_consensus(as.cl_ensemble(ensemble), method = "SE", control = list(nruns = 50))
    ensemble.labs <- c(cl_class_ids(cons))
    Ensemble <- c(Ensemble, adjustedRandIndex(trueLabels, ensemble.labs))
    Multiview <- c(Multiview, adjustedRandIndex(trueLabels, res@Labels$All))
    EnsembleNumClust <- c(EnsembleNumClust, length(unique(ensemble.labs)))
    MultiviewNumClust <- c(MultiviewNumClust, length(unique(res@Labels$All)))
  }
}
dats <- c("5000 Cells, 4 Clusters", "5000 Cells, 8 Clusters", "5000 Cells, 16 Clusters", 
          "10000 Cells, 4 Clusters", "10000 Cells, 8 Clusters", "10000 Cells, 16 Clusters",
          "15000 Cells, 4 Clusters", "15000 Cells, 8 Clusters", "15000 Cells, 16 Clusters")
result_table_balanced <- data.frame(Data = dats, Ensemble = Ensemble, Multiview = Multiview, EnsembleNumClust = EnsembleNumClust, 
                                    MultiviewNumClust = MultiviewNumClust)
write_csv(result_table_balanced, '/home/elijah/Documents/Clusterability/Results/Thesis/Simulation/Setup_1/Significance_Balanced.csv')

# unbalanced
Ensemble <- c()
Multiview <- c()
EnsembleNumClust <- c()
MultiviewNumClust <- c()
for(i in seq(5000, 15000, 5000)){
  for(j in c(4, 8, 16)){
    message(paste("Now running for numCells:", i))
    message(paste("Now running for numClusters:", j))
    sce <- Setup_1(numCells = i, numClusters = j, size = "unbalanced")
    meta.data <- colData(sce)
    trueLabels <- as.numeric(as.factor(meta.data$Group))
    
    # Run singnificance testing
    res <- RunSignificanceTest(counts(sce))
    
    # compute results
    ensemble.mat <- cbind(res@Labels$Seurat, res@Labels$Scrna, res@Labels$SCTransform)
    ensemble <- apply(ensemble.mat, MARGIN = 2, as.cl_partition)
    cons <- cl_consensus(as.cl_ensemble(ensemble), method = "SE", control = list(nruns = 50))
    ensemble.labs <- c(cl_class_ids(cons))
    Ensemble <- c(Ensemble, adjustedRandIndex(trueLabels, ensemble.labs))
    Multiview <- c(Multiview, adjustedRandIndex(trueLabels, res@Labels$All))
    EnsembleNumClust <- c(EnsembleNumClust, length(unique(ensemble.labs)))
    MultiviewNumClust <- c(MultiviewNumClust, length(unique(res@Labels$All)))
  }
}
dats <- c("5000 Cells, 4 Clusters", "5000 Cells, 8 Clusters", "5000 Cells, 16 Clusters", 
          "10000 Cells, 4 Clusters", "10000 Cells, 8 Clusters", "10000 Cells, 16 Clusters",
          "15000 Cells, 4 Clusters", "15000 Cells, 8 Clusters", "15000 Cells, 16 Clusters")
result_table_unbalanced <- data.frame(Data = dats, Ensemble = Ensemble, Multiview = Multiview, EnsembleNumClust = EnsembleNumClust, 
                                    MultiviewNumClust = MultiviewNumClust)
write_csv(result_table_unbalanced, '/home/elijah/Documents/Clusterability/Results/Thesis/Simulation/Setup_1/Significance_Unbalanced.csv')


# Simulation setup 2
Ensemble <- c()
Multiview <- c()
EnsembleNumClust <- c()
MultiviewNumClust <- c()
for(i in seq(0, 0.5, 0.01)){
    message(paste("Now running for separability value:", i))
    sce <- Setup_2(separability = i)
    meta.data <- colData(sce)
    trueLabels <- as.numeric(as.factor(meta.data$Group))
    
    # Run singnificance testing
    res <- RunSignificanceTest(counts(sce))
    
    # compute results
    ensemble.mat <- cbind(res@Labels$Seurat, res@Labels$Scrna, res@Labels$SCTransform)
    ensemble <- apply(ensemble.mat, MARGIN = 2, as.cl_partition)
    cons <- cl_consensus(as.cl_ensemble(ensemble), method = "SE", control = list(nruns = 50))
    ensemble.labs <- c(cl_class_ids(cons))
    Ensemble <- c(Ensemble, adjustedRandIndex(trueLabels, ensemble.labs))
    Multiview <- c(Multiview, adjustedRandIndex(trueLabels, res@Labels$All))
    EnsembleNumClust <- c(EnsembleNumClust, length(unique(ensemble.labs)))
    MultiviewNumClust <- c(MultiviewNumClust, length(unique(res@Labels$All)))
}
result_table <- data.frame(Separability = seq(0, 0.5, 0.01), Ensemble = Ensemble, Multiview = Multiview, 
                                    EnsembleNumClust = EnsembleNumClust, MultiviewNumClust = MultiviewNumClust)
write_csv(result_table, '/home/elijah/Documents/Clusterability/Results/Thesis/Simulation/Setup_2/Significance_Separability.csv')
# plot the results
d.ari <- dplyr::select(result_table, c("Separability", "Ensemble", "Multiview"))
d.ari.nelted <- melt(d.ari, id.vars = c("Separability"), variable.name = "Method", value.name = "ARI")
d.numclust <- dplyr::select(result_table, c("Separability", "EnsembleNumClust", "MultiviewNumClust"))
colnames(d.numclust) <- c("Separability", "Ensemble", "Multiview")
d.numclust.melted <-  melt(d.numclust, id.vars = c("Separability"), variable.name = "Method", value.name = "numClusters")

p.ari <- ggplot(d.ari.nelted, aes(Separability, ARI)) + geom_point(aes(color = Method))+
  ggtitle("ARI as a function of Separability") + facet_grid(. ~ Method) +
  theme(plot.title = element_text(hjust = 0.5, size=15), axis.title=element_text(size=14)) 

p.numclust <- ggplot(d.numclust.melted, aes(Separability, numClusters)) + geom_point(aes(color = Method))+
  ggtitle("Number of Clusters as a function of Separability") + facet_grid(. ~ Method) +
  theme(plot.title = element_text(hjust = 0.5, size=15), axis.title=element_text(size=14)) 


# Simulation setup 3
Ensemble <- c()
Multiview <- c()
EnsembleNumClust <- c()
MultiviewNumClust <- c()
for(i in seq(0, 5, 1)){
    message(paste("Now running for dropout rate:", i))
    sce <- Setup_3(dropout.mid = i)
    meta.data <- colData(sce)
    trueLabels <- as.numeric(as.factor(meta.data$Group))
    
    # Run singnificance testing
    res <- RunSignificanceTest(counts(sce))
    
    # compute results
    ensemble.mat <- cbind(res@Labels$Seurat, res@Labels$Scrna, res@Labels$SCTransform)
    ensemble <- apply(ensemble.mat, MARGIN = 2, as.cl_partition)
    cons <- cl_consensus(as.cl_ensemble(ensemble), method = "SE", control = list(nruns = 50))
    ensemble.labs <- c(cl_class_ids(cons))
    Ensemble <- c(Ensemble, adjustedRandIndex(trueLabels, ensemble.labs))
    Multiview <- c(Multiview, adjustedRandIndex(trueLabels, res@Labels$All))
    EnsembleNumClust <- c(EnsembleNumClust, length(unique(ensemble.labs)))
    MultiviewNumClust <- c(MultiviewNumClust, length(unique(res@Labels$All)))
}
result_table <- data.frame(DropoutRate = seq(0, 5, 1), Ensemble = Ensemble, Multiview = Multiview, 
                                    EnsembleNumClust = EnsembleNumClust, MultiviewNumClust = MultiviewNumClust)
write_csv(result_table, '/home/elijah/Documents/Clusterability/Results/Thesis/Simulation/Setup_3/Significance_Sparsity.csv')





















