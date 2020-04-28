# P <- mean(D <= replicate(B, dip(runif(n))))

a <- 0
b <- 1
d <- 100
n <- 1000
# samps <- cov(t(replicate(n, runif(d,a,b))) # draw samples
# d <- cov(t(samps))  # get the sample covariance matrix

# t <- replicate(1000, dip.test(parallelDist((replicate(n, runif(d,a,b)))))$p.value)

D <- as.data.frame(replicate(100000, dip(runif(10000))))
colnames(D) <- c("Dip")
p.dip <- ggplot(D, aes(x=Dip)) + geom_histogram() + 
  labs(title = "Dip Stat in for X in U[0,1]") + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90))
D.pvalue <- as.data.frame(replicate(100000, dip.test(runif(10000))$p.value)) 
colnames(D.pvalue) <- c("Pvalue")
p.pval <- ggplot(D.pvalue, aes(x=Pvalue)) + geom_histogram() + 
  labs(title = "Dip Stat P-value for X in U[0,1]") + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90))

p.final <- plot_grid(p.dip, p.pval, labels = c("A", "B"))
