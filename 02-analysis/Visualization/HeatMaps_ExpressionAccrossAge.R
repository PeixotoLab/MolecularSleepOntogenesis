
# CM, Oct 2022, plotting of the fold changes after SD in selected genes from functional cluster
# enrichment (DAVID)

library(grid)
library(dplyr)
library(plyr)
library(RColorBrewer)
library(forcats)
library(ggplot2)  
library(shades)
library(gridExtra)
library(gplots)
library(pheatmap)
library(here)

#  first load the normalized count matrix.
fileP16 <-here("Data","RUVOutput","HCvsSD_P16_ALL.txt")
fileP24 <-here("Data","RUVOutput","HCvsSD_P24_ALL.txt")
fileP30 <-here("Data","RUVOutput","HCvsSD_P30_ALL.txt")
fileP90 <-here("Data","GRUVOutput","HCvsSD_P90_ALL.txt")

# read the files in and change the colnames to make sure you do not mix them up
P16 <- read.table(fileP16, row.names = 1, header = TRUE)
colnames(P16) <- c("logFC_P16","logCPM_P16","F_P16","PValue_P16","FDR_P16")   
P16$ID <-rownames(P16)

P24 <- read.table(fileP24, row.names = 1, header = TRUE)
colnames(P24) <- c("logFC_P24","logCPM_P24","F_P24","PValue_P24","FDR_P24")   
P24$ID <-rownames(P24)

P30 <- read.table(fileP30, row.names = 1, header = TRUE)
colnames(P30) <- c("logFC_P30","logCPM_P30","F_P30","PValue_P30","FDR_P30")   
P30$ID <-rownames(P30)

P90 <- read.table(fileP90, row.names = 1, header = TRUE)
colnames(P90) <- c("logFC_P90","logCPM_P90","F_P90","PValue_P90","FDR_P90")  
P90$ID <-rownames(P90)
# make one large df from all ages 

masterdf <- join_all(list(P16,P24,P30,P90), by = 'ID', type = 'full')
row.names(masterdf) <-masterdf$ID
masterdf$ID <- NULL

#  define color range, make blue red and white
colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")) )(255)

# load the gene list of interest
#  First, genes with foldchanges in the Upper Quartile in Function Terms Up and Down, in the P24/P30/P90 group (3 for each direction)
# note: UP list with genes from 6 terms, DOWN list with genes from 2 terms
file1 <-here("Data","Genes4HeatMap","P24_P30_P90_UP_UQ.txt")
P24_P30_P90_UP <- read.table(file1,header = TRUE)

file2 <-here("Data","Genes4HeatMap","P24_P30_P90_DOWN_UQ.txt")
P24_P30_P90_DOWN <- read.table(file2,header = TRUE)

# alternative approach use only genes with terms that are at least 4 fold enriched
# note: UP list with genes from 2 terms, DOWN list with genes from 2 terms
# file1 <-here("Data","Genes4HeatMap","P24_P30_P90_UP_4FoldChange.txt")
# P24_P30_P90_UP <- read.table(file1,header = TRUE)
# 
# file2 <-here("Data","Genes4HeatMap","P24_P30_P90_DOWN_4FoldChange.txt")
# P24_P30_P90_DOWN <- read.table(file2,header = TRUE)



P24_P30_P90 <- rbind(P24_P30_P90_UP,P24_P30_P90_DOWN)
# make a new dataframe and find the foldChanges for each age group 

index1<- match(P24_P30_P90[,1],rownames(masterdf))
dfnew1 <- masterdf[index1,]

#  plot the log2FC for all ages
#  first extract only those values and make sure to keep the rownames with the unique IDs
FCmat1 <- rownames(dfnew1)
FCmat1 <- cbind(dfnew1$logFC_P16,dfnew1$logFC_P24,dfnew1$logFC_P30,dfnew1$logFC_P90)
rownames(FCmat1) <- rownames(dfnew1)
#  add the gene symbol in for plotting
index2 <-match(P24_P30_P90[,1],rownames(FCmat1))
symbol <-P24_P30_P90[index2,2]
rownames(FCmat1)<-symbol
colnames(FCmat1)<-c("P16","P24","P30","P90")
FCmat1 <- FCmat1[,2:4]
FCmat1 <- unique(FCmat1)

# plot as heatmap, euclidean distance with Ward's linkage
# euclidean method: \sqrt{\sum (i)\lef(x(i)-y(i)\right2}
# canberra method: Σ(|ai – bi|/(|ai|+|bi|)) absolute distance between two vectors divided by the sum of absolute of both vectors
# manhattan method: Σ|ai – bi| absolute distance between two vectors
heatmap.2(FCmat1, col = colors, Rowv=TRUE,Colv=NULL,
          distfun = function(x) dist(x, method="euclidean"),
          hclustfun = function(x) hclust(x, method="ward.D2"),
          dendrogram="row",scale="none", trace=c("none"), 
          key=TRUE,density.info=c("none"), key.title= NA, key.xlab="Log2 Foldchange SD/HC",main="DEG P24/P30/P90")

# pheatmap(FCmat1,color=colors,kmeans_k=NA, scale="column",cluster_rows=TRUE, cluster_cols=FALSE)


# # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # 

#  Group1, All Ages, genes with foldchanges in the Upper Quartile in Function Terms Up and Down, in the P24/P30/P90 group (3 for each direction)
# note: 1 term for the UP list, 3 terms for the DOWN list
file1 <-here("Data","Genes4HeatMap","AllAges_UP_UQ.txt")
AllAges_UP <- read.table(file1,header = TRUE)

file2 <-here("Data","Genes4HeatMap","AllAges_DOWN_UQ.txt")
AllAges_DOWN <- read.table(file2,header = TRUE)

# alternative approach use only genes with terms that are at least 4 fold enriched
# # note, the UP file does not change, the DOWN files comes from 2 terms
# file1 <-here("Data","Genes4HeatMap","AllAges_UP_4FoldChange.txt")
# AllAges_UP <- read.table(file1,header = TRUE)
# 
# file2 <-here("Data","Genes4HeatMap","AllAges_DOWN_4FoldChange.txt")
# AllAges_DOWN <- read.table(file2,header = TRUE)

AllAges <- rbind(AllAges_UP,AllAges_DOWN)
# make a new dataframe and find the foldChanges for each age group 

index1<- match(AllAges[,1],rownames(masterdf))
dfnew4 <- masterdf[index1,]

#  plot the log2FC for all ages
#  first extract only those values and make sure to keep the rownames with the unique IDs
FCmat4 <- rownames(dfnew4)
FCmat4 <- cbind(dfnew4$logFC_P16,dfnew4$logFC_P24,dfnew4$logFC_P30,dfnew4$logFC_P90)
rownames(FCmat4) <- rownames(dfnew4)
#  add the gene symbol in for plotting
index2 <-match(AllAges[,1],rownames(FCmat4))
symbol <-AllAges[index2,2]
rownames(FCmat4)<-symbol
colnames(FCmat4)<-c("P16","P24","P30","P90")

FCmat4 <- unique(FCmat4)
# plot as heatmap

# euclidean method: \sqrt{\sum (i)\lef(x(i)-y(i)\right2}
# canberra method: Σ(|ai – bi|/(|ai|+|bi|)) absolute distance between two vectors divided by the sum of absolute of both vectors
# manhattan method: Σ|ai – bi| absolute distance between two vectors
#  I like the canberra distance method the most
heatmap.2(FCmat4, col = colors, Rowv=TRUE,Colv=NULL, 
          distfun = function(x) dist(x, method="euclidean"),
          hclustfun = function(x) hclust(x, method="ward.D2"),
          dendrogram="row",scale="none", trace=c("none"), 
          key=TRUE,density.info=c("none"), key.title= NA, key.xlab="Log2 Foldchange SD/HC",main="DEG ALL ages ")


#  k-means clustering, using a k=3, because there genes come from around 3 therms mostly
kmeansObj <- kmeans(FCmat4, centers = 3)


heatmap(FCmat4 )
# =======================================
sink("sessionInfo_heatmaps.txt")
sessionInfo()
sink()
