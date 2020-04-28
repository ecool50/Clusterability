# res <- lapply(as.data.frame(dat), function(x){ length(which(x==0))/length(x)})
# Zero_Fraction <- as.numeric(unlist(res))

# dropout.mid = 0
sce <- Setup_3()
obj <- CreateModalityObject(counts(sce))
obj <- PreprocessObject(obj)
dat <- counts(sce)
res <- lapply(as.data.frame(dat), function(x){ length(which(x==0))/length(x)})
Zero_Fraction.0 <- as.numeric(unlist(res))
x <- as.data.frame(obj@Seurat$Sig.PCs.Data)
p.0 <- ggplot(x, aes(V1, V2, color = Zero_Fraction.0)) + geom_point() + ggtitle("Low Sparsity") + xlab('PC1') + ylab("PC2") +
  theme(plot.title = element_text(hjust = 0.5, size=15), axis.title=element_text(size=14)) + labs(color = "Zero Fraction")
  

# dropout.mid = 1
sce <- Setup_3(dropout.mid = 1)
obj <- CreateModalityObject(counts(sce))
obj <- PreprocessObject(obj)
dat <- counts(sce)
res <- lapply(as.data.frame(dat), function(x){ length(which(x==0))/length(x)})
Zero_Fraction.1 <- as.numeric(unlist(res))
x <- as.data.frame(obj@Seurat$Sig.PCs.Data)
p.1 <- ggplot(x, aes(V1, V2, color = Zero_Fraction.1)) + geom_point() + ggtitle("dropout.mid = 1") + xlab('PC1') + ylab("PC2") +
theme(plot.title = element_text(hjust = 0.5, size=15), axis.title=element_text(size=14))  + labs(color = "Zero Fraction")


# dropout.mid = 3
sce <- Setup_3(dropout.mid = 3)
obj <- CreateModalityObject(counts(sce))
obj <- PreprocessObject(obj)
dat <- counts(sce)
res <- lapply(as.data.frame(dat), function(x){ length(which(x==0))/length(x)})
Zero_Fraction.3 <- as.numeric(unlist(res))
x <- as.data.frame(obj@Seurat$Sig.PCs.Data)
p.3 <- ggplot(x, aes(V1, V2, color = Zero_Fraction.3)) + geom_point() + ggtitle("dropout.mid = 3") + xlab('PC1') + ylab("PC2") +
theme(plot.title = element_text(hjust = 0.5, size=15), axis.title=element_text(size=14)) + labs(color = "Zero Fraction")

# dropout.mid = 5
sce <- Setup_3(dropout.mid = 5)
obj <- CreateModalityObject(counts(sce))
obj <- PreprocessObject(obj)
dat <- counts(sce)
res <- lapply(as.data.frame(dat), function(x){ length(which(x==0))/length(x)})
Zero_Fraction.5 <- as.numeric(unlist(res))
x <- as.data.frame(obj@Seurat$Sig.PCs.Data)
p.5 <- ggplot(x, aes(V1, V2, color = Zero_Fraction.5)) + geom_point() + ggtitle("High Sparsity") + xlab('PC1') + ylab("PC2") +
theme(plot.title = element_text(hjust = 0.5, size=15), axis.title=element_text(size=14)) + labs(color = "Zero Fraction")

p.final <- plot_grid(p.0, p.1, p.3, p.5, labels = c("A", "B", "C", "D"))
