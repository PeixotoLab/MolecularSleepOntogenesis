#  script to run a hierarchial k-means clustering for the developmental data, AFTER RUVs.
#  Christine Muheim, June 2022

# Load dependencies
library("ggplot2")
library("ggdendro")
library("tidyr")
library("grid")
library("RColorBrewer")
library("gplots")
library("dendextend")
library("here")
# Read in data
file1 <-here("Data","GeneListALLages_Filtered.txt")
dataALL <- read.table(file1, row.names = 1, header = TRUE)
#  this data is filtered, and all the genes, not yet after differential expression!!
# define the groups, might help to lable later on
x <- as.factor(rep(c("HC16","HC24","HC30","HC90","SD16","SD24","SD30","SD90"), c(5,8,5,8,5,8,5,8))) 
names(x) <- colnames(dataALL)

# plot variance vs mean to see how variable the data is, I used this tutorial to guide me through this all: 
# https://2-bitbio.com/2017/10/clustering-rnaseq-data-using-k-means.html

z <- dataALL
z_var <- apply(z, 1, var)
z_mean <- apply(z, 1, mean)
plot(log2(z_mean), log2(z_var), pch='.')
#  plot the varibility of the data, just to see
abline(h=log2(175), col='red')
abline(v=log2(175), col='red')
text(x=13,y=23, labels="variance >75 &\n mean > 75", col='red')


scaledata <- t(scale(t(z))) # Centers and scales data. This basically does this: $ data(scaled) = data(i)-mean(data) / sd(data) $, 
# meaning that the distribution has mean 0 and standard deviation (sd) 1.

#  cluster samples to identify out layers

hc <- hclust(as.dist(1-cor(scaledata, method="spearman")), method="average") # Clusters columns by Spearman correlation.
TreeC = as.dendrogram(hc, method="average")
plot(TreeC,
     main = "Sample Clustering",
     ylab = "Height")

# heatmap.2(as.matrix(dataALL[,-1], trace = "none", na.color = "Green"))

# =======================================
# =======================================
# =======================================
sink("sessionInfo_Clustering.txt")
sessionInfo()
sink()


small_iris <- iris[c(1, 51, 101, 2, 52, 102), ]
dend <- as.dendrogram(hclust(dist(small_iris[,-5])))
labels_colors(dend)
colors_to_use <- as.numeric(small_iris[,5])
