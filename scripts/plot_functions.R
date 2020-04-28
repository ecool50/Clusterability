library(tidyverse)
library(gridExtra)
setwd("/home/elijah/Documents/Clusterability")

#read in the files for plotting
data.1 <- read.csv("GF_Frequencies.csv", header = T)
data.2 <- read.csv("Citro_Frequencies.csv", header = T)
data.3<- read.csv("BR_Frequencies.csv", header = T)

#plot the data
p1 <- ggplot(data.1, aes(x=Frequency)) + 
  geom_histogram(colour="black", fill="lightblue") + ggtitle("GF") + geom_vline(aes(xintercept=mean(Frequency)),
                                                                      color="black", linetype="dashed", size=1)

p2 <- ggplot(data.2, aes(x=Frequency)) + 
  geom_histogram(colour="black", fill="lightblue") + ggtitle("Citro") + geom_vline(aes(xintercept=mean(Frequency)),
                                                                           color="black", linetype="dashed", size=1)

p3 <- ggplot(data.3, aes(x=Frequency)) + 
  geom_histogram(colour="black", fill="lightblue") + ggtitle("BR") + geom_vline(aes(xintercept=mean(Frequency)),
                                                                      color="black", linetype="dashed", size=1)

#plot them all on one grid
p.final <- grid.arrange(p1, p2, p3, ncol=3)
