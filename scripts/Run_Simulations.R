# Run simulation setup 1
# For the balanced clusters
Log <- c()
NegBinom <- c()
Multinom <- c()
for(i in seq(5000, 15000, 5000)){
  for(j in c(4, 8, 16)){
    message(paste("Now running for numCells:", i))
    message(paste("Now running for numClusters:", j))
    sce <- Setup_1(numCells = i, numClusters = j)
    res <- t(as.data.frame(RunModalityTest(counts(sce))))
    Log <- c(Log, res[[1]])
    Multinom <- c(Multinom, res[[2]])
    NegBinom <- c(NegBinom, res[[3]])
  }
}
dats <- c("5000 Cells, 4 Clusters", "5000 Cells, 8 Clusters", "5000 Cells, 16 Clusters", 
          "10000 Cells, 4 Clusters", "10000 Cells, 8 Clusters", "10000 Cells, 16 Clusters",
          "15000 Cells, 4 Clusters", "15000 Cells, 8 Clusters", "15000 Cells, 16 Clusters")
result_table_balanced <- data.frame(Data = dats, Log = Log, Multinom = Multinom, NegBinom = NegBinom)
write_csv(result_table_balanced, '/home/elijah/Documents/Clusterability/Results/Thesis/Simulation/Setup_1/Balanced.csv')

# for the inbalanced Clusters
Log <- c()
NegBinom <- c()
Multinom <- c()
for(i in seq(5000, 15000, 5000)){
  for(j in c(4, 8, 16)){
    message(paste("Now running for numCells:", i))
    message(paste("Now running for numClusters:", j))
    sce <- Setup_1(numCells = i, numClusters = j, size = "unbalanced")
    res <- t(as.data.frame(RunModalityTest(counts(sce))))
    Log <- c(Log, res[[1]])
    Multinom <- c(Multinom, res[[2]])
    NegBinom <- c(NegBinom, res[[3]])
  }
}
dats <- c("5000 Cells, 4 Clusters", "5000 Cells, 8 Clusters", "5000 Cells, 16 Clusters", 
          "10000 Cells, 4 Clusters", "10000 Cells, 8 Clusters", "10000 Cells, 16 Clusters",
          "15000 Cells, 4 Clusters", "15000 Cells, 8 Clusters", "15000 Cells, 16 Clusters")
result_table_unbalanced <- data.frame(Data = dats, Log = Log, Multinom = Multinom, NegBinom = NegBinom)
write_csv(result_table_unbalanced, '/home/elijah/Documents/Clusterability/Results/Thesis/Simulation/Setup_1/Unbalanced.csv')

# Run simulation setup 2
Log <- c()
NegBinom <- c()
Multinom <- c()
for (i in seq(0, 0.5, 0.01)){
  message(paste("Now running for separability value:", i))
  sce <- Setup_2(separability = i)
  res <- t(as.data.frame(RunModalityTest(counts(sce))))
  Log <- c(Log, res[[1]])
  Multinom <- c(Multinom, res[[2]])
  NegBinom <- c(NegBinom, res[[3]])
}

# plot the results
result_table <- data.frame(Separability = seq(0, 0.5, 0.01), Log = Log, Multinom = Multinom, NegBinom = NegBinom)
dip_vals <- melt(result_table, id.vars = "Separability")
colnames(dip_vals) <- c("Separability", "Model", "value")
dip_plot <- ggplot2::ggplot(dip_vals) + geom_point(aes(Separability, value, color = Model)) +
  xlab("Cluster Separability") + ylab("Dip Pvalue") + theme_gray() +
  ggtitle("Clusterability as a Function of Cluster Seperability") +
  facet_grid(. ~ Model) +
  geom_hline(yintercept=0.05, linetype="dashed", color = "red") + 
  theme(plot.title = element_text(hjust = 0.5, size=15), axis.title=element_text(size=14))

write_csv(result_table, '/home/elijah/Documents/Clusterability/Results/Thesis/Simulation/Setup_2/Separability.csv')
pdf("/home/elijah/Documents/Clusterability/plots/Thesis_Plots/Simulations/Setup_2/Separability.pdf", width = 12, height = 8)
plot(dip_plot)
dev.off()

# Run simulation setup 3
Log <- c()
NegBinom <- c()
Multinom <- c()
for (i in seq(0, 5, 1)){
  message(paste("Now running for dropout value:", i))
  sce <- Setup_3(dropout.mid = i)
  res <- t(as.data.frame(RunModalityTest(counts(sce))))
  Log <- c(Log, res[[1]])
  Multinom <- c(Multinom, res[[2]])
  NegBinom <- c(NegBinom, res[[3]])
}

result_table <- data.frame(DropoutRate = seq(0, 5, 1), Sparsity = c(0.25, 0.34, 0.46, 0.61, 0.75, 0.85), Log = Log, 
                           Multinom = Multinom, NegBinom = NegBinom)
# plot the results
p.sparsity <- dat.melted %>% ggplot(aes(x = DropoutRate, y = value)) + geom_point(aes(color = Sparsity)) + 
  facet_grid(. ~ variable) + geom_hline(yintercept=0.05, linetype="dashed", color = "red") + 
  xlab("dropout.mid") + ylab("Dip Pvalue") + theme_gray() +
  ggtitle("Clusterability as a Function of Data Sparsity") +
  theme(plot.title = element_text(hjust = 0.5, size=15), axis.title=element_text(size=14))
pdf("/home/elijah/Documents/Clusterability/plots/Thesis_Plots/Simulations/Setup_3/Sparsity.pdf", width = 12, height = 8)
plot(p.sparsity)
dev.off()
write_csv(result_table, '/home/elijah/Documents/Clusterability/Results/Thesis/Simulation/Setup_3/Sparsity.csv')
  
