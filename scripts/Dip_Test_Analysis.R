library(diptest)
library(ggplot2)


# CONSTANT PARAMETERS
sig1 <- log(3)
sig2 <- log(3)
cpct <- 0.5
N=100

#CREATING BIMOD DISTRIBUTION
bimodalDistFunc <- function (n,cpct, mu1, mu2, sig1, sig2) {
  y0 <- rlnorm(n,mean=mu1, sd = sig1)
  y1 <- rlnorm(n,mean=mu2, sd = sig2)
  
  flag <- rbinom(n,size=1,prob=cpct)
  y <- y0*(1 - flag) + y1*flag 
}

#DIP TEST
DIP_TEST <- function(bimodalData) {
  TEST <- dip.test(bimodalData, simulate.p.value = TRUE)
    # modetest(bimodalData, method = "CH", B = 100)
    # dip.test(bimodalData)
  return(list(TEST$statistic[[1]], TEST$p.value[[1]]))   # return(TEST$p.value[[1]])    to get the p value
}
DIP_TEST(bimodalData)


# SIMULATION
exp_mu1 = 1
max_exp_mu2 = 100
intervStep = 100
repPerInt = 10

# single distibutions
expMu2Value <- c()
bimodalData <- c()
mu1 <- log(exp_mu1)   
mu2 <- log(exp_mu1)
bimodalData <- c(bimodalData,log(bimodalDistFunc(n=N,cpct,mu1,mu2, sig1,sig2)))
expMu2Value <- c(expMu2Value,rep(exp_mu1,length(log(bimodalDistFunc(n=N,cpct,mu1,mu2, sig1,sig2)))))

mu1 <- log(exp_mu1)   
mu2 <- log(max_exp_mu2)
bimodalData <- c(bimodalData,log(bimodalDistFunc(n=N,cpct,mu1,mu2, sig1,sig2)))
expMu2Value <- c(expMu2Value,rep(max_exp_mu2,length(log(bimodalDistFunc(n=N,cpct,mu1,mu2, sig1,sig2)))))

mu1 <- log(exp_mu1)   
mu2 <- log(trunc((max_exp_mu2-exp_mu1)/2+1))
bimodalData <- c(bimodalData,log(bimodalDistFunc(n=N,cpct,mu1,mu2, sig1,sig2)))
expMu2Value <- c(expMu2Value,rep(trunc((max_exp_mu2-exp_mu1)/2+1),length(log(bimodalDistFunc(n=N,cpct,mu1,mu2, sig1,sig2)))))

tableExamples <- data.frame(expMu2Value,bimodalData)
tableExamples$expMu2Value <- as.factor(tableExamples$expMu2Value)
ExamplePlot <- ggplot(tableExamples)+
  geom_histogram(aes(bimodalData),color='white')+
  ylab("Count")+
  xlab("")+
  facet_wrap(~expMu2Value)+
  ggtitle("Intensity of bimodularity")

# calculation of the dip test index
exp_mu2Int = seq(from=exp_mu1,to=max_exp_mu2,length.out=intervStep)
expmu2Vec = c()
dipStat = c()
dipPval = c()
testDone = c()
for(exp_mu2 in exp_mu2Int){
  mu1 <- log(exp_mu1)   
  mu2 <- log(exp_mu2)
  for(rep in 1:repPerInt){
    bimodalData <- log(bimodalDistFunc(n=N,cpct,mu1,mu2, sig1,sig2))
    diptestone = DIP_TEST(bimodalData)[[1]]
    dip_pval = DIP_TEST(bimodalData)[[2]]
    expmu2Vec = c(expmu2Vec,exp_mu2)
    dipStat = c(dipStat,diptestone)
    dipPval = c(dipPval, dip_pval)
    testDone = c(testDone,"diptest")
  }
}
table = data.frame(expmu2Vec,dipStat,testDone)

IndexPlot <- ggplot(table)+
  geom_point(aes(expmu2Vec,dipStat,color=testDone))+
  ylab("Index")+
  xlab("Intensity of Bimodularity")+
  scale_color_discrete(name="Test")

pval_table <- data.frame(expmu2Vec, dipPval, testDone)

Pval_Plot = ggplot(table)+
  geom_point(aes(expmu2Vec,dipPval,color=testDone))+
  ylab("P value")+
  xlab("Intensity of Bimodularity")+
  scale_color_discrete(name="Test")
# pdf("/home/elijah/Documents/Clusterability/Results/Dip_Test_Analysis.pdf")
# ExamplePlot
# IndexPlot
# Pval_Plot
# dev.off()