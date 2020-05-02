res <- clusterlab(centers = 3, features = 3, sdvec = c(1,1,1), numbervec = c(300,300, 300))
x <- t(as.matrix(res$synthetic_data))
hc <- hclust(parallelDist(x), method = "ward.D2")

dend <- as.dendrogram(hc)
# dend <- dend %>% set("branches_k_color", k = 2)
# dend <- dend %>% set("hang",0.1)

ggd1 <- as.ggdend(dend)

n <- nrow(x)
idx_hc <- .idx_hc(hc, n)

idx_vals <- idx_hc[2, ]
x1 <- x[idx_vals[[1]], ]
x2 <- x[idx_vals[[2]], ]
x_comb <- rbind(x1,x2)[, 1:2]
assignments <- c(rep(1,nrow(x1)), rep(2, nrow(x2)))

labs <- rownames(x_comb)

dend %>% set("by_labels_branches_col", value = labs) %>% 
  plot(main = "node", leaflab = "none")

# prun the tree

dend %>% set("labels_colors") %>%
  prune(labs) %>% 
  plot(main = "Prunned tree")
# proj <- fpc::discrcoord(xd = x_comb, clvecd = assignments)

res.fisher <- ComputeFisher(x_comb = x_comb, assignments = assignments)


# plot the results
x <- as.data.frame(x_comb)
labels <- as.factor(assignments)
Clusters <- as.factor(cutree(dend, 2))

# plot the pca
p.pca <- as.data.frame(x_comb) %>% ggplot(aes(V1, V2)) + geom_point(aes(color = Clusters)) + ggtitle("2 Simulated Clusters") +
  theme(plot.title = element_text(hjust = 0.5, size=15), axis.title=element_text(size=14))

# plot the dendrogram
p.dend <- ggplot(ggd1, labels = F) + theme_gray() + ggtitle("Ward's Clustering")  + xlab("") + ylab("Linkage") +
  theme(plot.title = element_text(hjust = 0.5, size=15), axis.title=element_text(size=14))

# plot fischer Discriminant line
p.fisher <- ggplot(x, aes(x=V1,y=V2)) + geom_point(aes(color = labels)) + geom_abline(aes(intercept = res.fisher$intercept, slope = res.fisher$slope)) +
  ggtitle("Fisher Discriminant") + theme(plot.title = element_text(hjust = 0.5, size=15), axis.title=element_text(size=14))

# plot the density
r <- density(res.fisher$proj)
df <- data.frame(x=r$x, y=r$y)
p.density <- df %>% ggplot(aes(x,y)) + geom_line() + xlab("Fischer Projection") + ylab("Density") + 
  ggtitle("Discriminant Projection Density(Dip pval = 0)") + theme(plot.title = element_text(hjust = 0.5, size=15), axis.title=element_text(size=14))

p.final <- plot_grid(p.pca, p.dend, p.fisher, p.density, labels = c("A", "B", "C", "D"))
  