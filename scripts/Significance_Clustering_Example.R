synthetic <- clusterlab(centers=2,centralcluster=TRUE, numbervec = c(25, 25))

x <- t(as.data.frame(synthetic$synthetic_data))

x.pca <- as.data.frame(irlba::prcomp_irlba(t(x), n = 2)$rotation)

dend <- as.matrix(x) %>%  scale %>% 
  parallelDist %>% hclust(method = "ward.D2") %>% as.dendrogram

dend <- dend %>% set("hang",0.1)

dend <- dend %>% set("nodes_pch", 19) %>%  # node point type
  set("nodes_cex", 2) %>%  # node point size
  set("nodes_col", "blue")

dend <- dend %>% set("leaves_pch", 19) %>%  # node point type
  set("leaves_cex", 2) %>%  # node point size
  set("leaves_col", "red") 

ggd1 <- as.ggdend(dend)

p.dend <- ggplot(ggd1, labels = F) + theme_gray() + ggtitle("Ward's Clustering")  + xlab("") + ylab("Linkage") +
  theme(plot.title = element_text(hjust = 0.5, size=15), axis.title=element_text(size=14))

Clusters <- as.factor(cutree(dend, 2))
p.pca <- x.pca %>% ggplot(aes(PC1, PC2)) + geom_point(aes(color = Clusters)) + ggtitle("2 Simulated Gaussian Clusters") +
  theme(plot.title = element_text(hjust = 0.5, size=15), axis.title=element_text(size=14))
p.final <- plot_grid(p.pca, p.dend, labels = c("A", "B"))

### Plot modality example
sce <- Setup_1()
obj <- CreateModalityObject(counts(sce), sparse = FALSE)
# add the seurat slot
obj <- PreprocessObject(obj,  nFeatures = 500)
# add the sctransform slot
obj <- PreprocessObject(obj,  method = "SCTransform", nFeatures = 500)
# add the scrna slot
obj <- PreprocessObject(obj, method = "Scrna", nFeatures = 500)

sce <- scater::normalize(sce)
p.pca <- plotPCA(sce, colour_by = "Group") + ggtitle(paste("Four Simulated Clusters")) + 
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90)) +  theme_grey() + 
  theme(plot.title = element_text(hjust = 0.5, size=15), axis.title=element_text(size=14))

# log
d1 <- coop::tcosine(obj@Seurat$Sig.PCs.Data)
d1 <- density((1 - d1)/2)
d1 <- data.frame(x = d1$x, y= d1$y)
p1 <- ggplot(d1, aes(x,y)) + geom_line(color = "blue") + theme_gray() + ggtitle("Log Normalisation") + 
  theme(plot.title = element_text(hjust = 0.5, size=15), axis.title=element_text(size=14)) +
  xlab("Cosine Distance") + ylab("Density")
  

# multinomial
d2 <- coop::tcosine(obj@Scrna$Sig.PCs.Data)
d2 <- density((1 - d2)/2)
d2 <- data.frame(x = d2$x, y= d2$y)
p2 <- ggplot(d2, aes(x,y)) + geom_line(color = "blue") + theme_gray() + ggtitle("Multinomial Normalisation") + 
  theme(plot.title = element_text(hjust = 0.5, size=15), axis.title=element_text(size=14)) +
  xlab("Cosine Distance") + ylab("Density")

# negBinom
d3 <- coop::tcosine(obj@SCTransform$Sig.PCs.Data)
d3 <- density((1 - d3)/2)
d3 <- data.frame(x = d3$x, y= d3$y)
p3 <- ggplot(d3, aes(x,y)) + geom_line(color = "blue") + theme_gray() + ggtitle("NegBinomial Normalisation") + 
  theme(plot.title = element_text(hjust = 0.5, size=15), axis.title=element_text(size=14)) +
  xlab("Cosine Distance") + ylab("Density")

p.final.modality <- plot_grid(p.pca, p1, p2, p3, labels = c("A", "B", "C", "D"))

