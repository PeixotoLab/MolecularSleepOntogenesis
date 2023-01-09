#  CM, March, 2022, import the DEG from the excel table and start getting the venn/euler diagrams together
# find the shared or unique gene sets of interest and export them (make sure to have EnsembleID, log2FC, FDR and pValue and Annotation in there too)


######### VENN diagrams
library(ggVennDiagram)
library(grid)
library(ggplot2)
library(readxl)
library(dplyr)
library(RColorBrewer)
library(eulerr)
library(here)
library(ggvenn)


sink("sessionInfo_Euler.txt")
sessionInfo()
sink()

####### use the exported data in the text files and only import the first column for this

####### load all the lists, these are the significantly changed genes, FDR 0.05
####### NOTE: these lists have been annotated and genes without annotation were removed, on average 2 genes per age group

outfile <- here("Data","GeneLists_RUVoutput","DEG_RUVk3_Significants_AllAges_May2022.xlsx")
P16 <- read_excel(outfile, sheet=1,
                  col_names = c("Gene_StableID","logFC_P16", "logCPM_P16", "F_P16", "PValue_P16", "FDR_P16","GeneName","GeneDescription"), skip =1)
P24 <- read_excel(outfile, sheet=2,
                  col_names = c("Gene_StableID","logFC_P24", "logCPM_P24", "F_P24", "PValue_P24", "FDR_P24","GeneName","GeneDescription"), skip =1)
P30 <- read_excel(outfile,sheet=3,
                  col_names = c("Gene_StableID", "logFC_P30", "logCPM_P30", "F_P30", "PValue_P30", "FDR_P30","GeneName","GeneDescription"), skip = 1)
P90 <- read_excel(outfile,sheet=4,
                  col_names = c("Gene_StableID", "logFC_P90", "logCPM_P90", "F_P90", "PValue_P90", "FDR_P90","GeneName","GeneDescription"), skip = 1)


# outfile2 <- "FinalDEGlistALL_Annotated"
# 
# P24_all <- read_excel(outfile, sheet=1,
#                       col_names = c("Gene_StableID","logFC_P24", "logCPM_P24", "F_P24", "PValue_P24", "FDR_P24"), skip =1)
# P30_all <- read_excel(outfile, sheet=2,
#                       col_names = c("Gene_StableID","logFC_P30", "logCPM_P30", "F_P30", "PValue_P30", "FDR_P30"), skip =1)
# P90_all <- read_excel(outfile, sheet=3,
#                       col_names = c("Gene_StableID","logFC_P90", "logCPM_P90", "F_P90", "PValue_P90", "FDR_P90"), skip =1)

colors <- brewer.pal(8, "Dark2")
colors2 <- brewer.pal(9,"Set1") 

color <- brewer.pal(4,"RdYlBu")
colorUP <- brewer.pal(4,"Reds")
colorDOWN <- brewer.pal(4,"Blues")

color1 <- adjust_transparency(color, alpha = 0.5)

####### use ENSEMBL ID to find what intersects
w <- P16$Gene_StableID
x <- P24$Gene_StableID
y <- P30$Gene_StableID
z <- P90$Gene_StableID

length(w)
length(x)
length(y)
length(z)

data <- list(A=w,B=x,C=y,D=z)

names(data) <- c("P16","P24","P30","P90")
ggvenn(data,
       fill_color = color,
       stroke_size = 0.25,
       set_name_size = 5,
       text_size =5,
       show_percentage = FALSE)

plot(euler(data,
           shape = "ellipse"),
           # shape = "circle") ,
           fills = color,
           alpha = 0.5,
           quantities = TRUE)



######### create a Venn/Euler diagram with the up and down regulated lists
P16_UP <- (P16[P16$logFC_P16>0,])
P16_DOWN <- (P16[P16$logFC_P16<0,])

P24_UP <- (P24[P24$logFC_P24>0,])
P24_DOWN <- (P24[P24$logFC_P24<0,])

P30_UP <- (P30[P30$logFC_P30>0,])
P30_DOWN <- (P30[P30$logFC_P30<0,])

P90_UP <- (P90[P90$logFC_P90>0,])
P90_DOWN <- (P90[P90$logFC_P90<0,])

w1 <- P16_UP$Gene_StableID
x1 <- P24_UP$Gene_StableID
y1 <- P30_UP$Gene_StableID
z1 <- P90_UP$Gene_StableID

w2 <- P16_DOWN$Gene_StableID
x2 <- P24_DOWN$Gene_StableID
y2 <- P30_DOWN$Gene_StableID
z2 <- P90_DOWN$Gene_StableID

dataUP <- list(A=w1,B=x1,C=y1,D=z1)
dataDOWN <- list(A=w2,B=x2,C=y2,D=z2)

names(dataUP) <- c("P16","P24","P30","P90")
ggvenn(dataUP,
      fill_color = colorUP,
      stroke_size = 0.25,
      set_name_size = 5,
      text_size =5,
      show_percentage = FALSE)

plot(euler(dataUP,
           shape = "ellipse"),
     fills = colorUP,
     alpha = 0.5,
     quantities = TRUE)



names(dataDOWN) <- c("P16","P24","P30","P90")
ggvenn(dataDOWN,
      fill_color = colorDOWN,
      stroke_size = 0.25,
      set_name_size = 5,
      text_size =5,
      show_percentage = FALSE)


plot(euler(dataDOWN,
           shape = "ellipse"),
     fills = colorDOWN,
     alpha = 0.5,
     quantities = TRUE)


# ################# get the lists of interest and export as txt file
# ###### first find all the genes in the class of interest
commonP24P30 <- intersect(P24$Gene_StableID,P30$Gene_StableID) ###shared between P24 and P30
commonP16P24P30 <- intersect(commonP24P30,P16$Gene_StableID) ###shared between P24 and P30 and P16, total of 5
commonALL <- intersect(commonP16P24P30,P90$Gene_StableID)###shared between all ages, total of 59
commonP24P30P90 <- intersect(commonP24P30,P90$Gene_StableID)###shared between P24, P30 and P90 total of 1563

commonP24P90 <- intersect(P24$Gene_StableID, P90$Gene_StableID) ###shared between P24 and P90
commonP30P90 <- intersect(P30$Gene_StableID, P90$Gene_StableID)###shared between P30 and P90

# all DEG in all the young ages
ALLyoung <- union(P24$Gene_StableID,P30$Gene_StableID)
ALLyoung <- union(P16$Gene_StableID,ALLyoung)
P24P90union <-union(P24$Gene_StableID,P90$Gene_StableID)
P30P90union <- union(P30$Gene_StableID,P90$Gene_StableID)
P30P90P16union <-union(P16$Gene_StableID,P30P90union)
P24P90P16union <-union(P16$Gene_StableID,P24P90union)
P24P30P90union <- union(P24$Gene_StableID,P30P90union)
#  get the groups that are unique to each age
P90only <- setdiff(P90$Gene_StableID,ALLyoung)
P30only <- setdiff(P30$Gene_StableID,P24P90P16union)
P24only <- setdiff(P24$Gene_StableID,P30P90P16union)
P16only <- setdiff(P16$Gene_StableID,P24P30P90union)

# shared between the P24 and P30 but not P90 
P24P30 <- setdiff(commonP24P30,P90$Gene_StableID)
P24P30only <- setdiff(P24P30,P16$Gene_StableID) # shared only between P24 and P30, not including the P16

#  shared between P24, P30 and P90 but NOT P16
exclusiveP24P30P90 <- setdiff(commonP24P30P90,P16$Gene_StableID)
# ###### now get the full list extracted, including FC, logCPM, F, P and FDR

indx1a <- match(commonALL,P16$Gene_StableID)
indx1b <- match(commonALL,P24$Gene_StableID)
indx1c <- match(commonALL,P30$Gene_StableID)
indx1d <- match(commonALL,P90$Gene_StableID)

SLEEPList <- cbind(P16[indx1a,],P24[indx1b,],P30[indx1c,],P90[indx1d,])# all genes that are DEG at all ages
file<-here("Data","OUTPUT_ListJune","DEG_AllAges.txt")
write.table(x = SLEEPList, file, sep = "\t")

indx2 <- match(P90only,P90$Gene_StableID)
OnlyP90 <-P90[indx2,] # all genes that are DEG only in P90
file<-here("Data","OUTPUT_ListJune","DEG_P90only.txt")
write.table(x = OnlyP90 , file, sep = "\t")

indx3 <- match(P30only,P30$Gene_StableID)
OnlyP30 <-P30[indx3,] # all genes that are DEG onlyin P30
file<-here("Data","OUTPUT_ListJune","DEG_P30only.txt")
write.table(x = OnlyP30 , file, sep = "\t")

indx4 <- match(P24only,P24$Gene_StableID)
OnlyP24 <-P24[indx4,] # all genes that are DEG onlyin P24
file<-here("Data","OUTPUT_ListJune","DEG_P24only.txt")
write.table(x = OnlyP24, file, sep = "\t")

indx5 <- match(P16only,P16$Gene_StableID)
OnlyP16 <-P16[indx5,] # all genes that are DEG onlyin P16
file<-here("Data","OUTPUT_ListJune","DEG_P16only.txt")
write.table(x = OnlyP16 , file, sep = "\t")

indxA <- match(exclusiveP24P30P90,P24$Gene_StableID) # all genes that are DEG in P24, P30 and P90 but not in P16
indxB <- match(exclusiveP24P30P90,P30$Gene_StableID) # all genes that are DEG in P24, P30 and P90 but not in P16
indxC <- match(exclusiveP24P30P90,P90$Gene_StableID) # all genes that are DEG in P24, P30 and P90 but not in P16

P24P30P90DEG_1 <-cbind(P24[indxA,],P30[indxB,],P90[indxC,])
P24P30P90DEG <- P24P30P90DEG_1[!duplicated(as.list(P24P30P90DEG_1))]

#  remove the duplicated lists late ron, there is a way in R but I have no time rigth now
file <- here("Data","OUTPUT_ListJune","DEG_P24P30P90.txt")
write.table(x = P24P30P90DEG , file, sep = "\t")

#  find all genes that are unique to P24 AND P30, not present in P16 or P90
indxI <- match(P24P30only,P24$Gene_StableID)
indxII <- match(P24P30only,P30$Gene_StableID)

P24_P30DEG <- cbind(P24[indxI,],P30[indxII,])
P24_P30DEG <- P24_P30DEG[!duplicated(as.list(P24_P30DEG))]
file <- here("Data","OUTPUT_ListJune","DEG_P24P30.txt")
write.table(x = P24_P30DEG , file, sep = "\t")

# =======================================
# =======================================
# =======================================
#   Euler diagram/ Gene lists for grant due June 5th, this is only to compare P30 and P90, do an all, up and down
# 
# dataUPgrant <- list(A=y1,B=z1)
# dataDOWNgrant <- list(A=y2,B=z2)
# dataGrant <-list(A=y,B=z)
# names(dataUPgrant) <- c("P30","P90")
# names(dataDOWNgrant) <- c("P30","P90")
# names(dataGrant) <- c("P30","P90")
# 
# plot(euler(dataGrant,
#            shape = "ellipse"),
#      fills = color,
#      alpha = 0.5,
#      quantities = TRUE
# )
# 
# plot(euler(dataDOWNgrant,
#            shape = "ellipse"),
#      fills = colorDOWN,
#      alpha = 0.5,
#      quantities = TRUE)
# 
# plot(euler(dataUPgrant,
#            shape = "ellipse"),
#      fills = colorUP,
#      alpha = 0.5,
#      quantities = TRUE)
# 
# # export only the shared genes for the P30/P90 interaction, also export the P30 onyl and the P90 only for this grant
# # find what P30 and P90 share
# commonP30P90 <- intersect(P30$Gene_StableID,P90$Gene_StableID)
# P30onlyGRN <- setdiff(P30$Gene_StableID,commonP30P90)
# P90onlyGRN <- setdiff(P90$Gene_StableID,commonP30P90)
# 
# indx1 <- match(P30onlyGRN,P30$Gene_StableID)
# indx2 <- match(P90onlyGRN,P90$Gene_StableID)
# indxA <- match(commonP30P90,P30$Gene_StableID) # all genes that are DEG in P24, P30 and P90 but not in P16
# indxB <- match(commonP30P90,P90$Gene_StableID)
# 
# OnlyP30Grn <- P30[indx1,]
# OnlyP90Grn <- P90[indx2,]
# 
# Shared_P30P90 <- cbind(P30[indxA,],P90[indxB,])
# duplicatedRows <- duplicated(colnames(Shared_P30P90))
# Shared_P30P90[!duplicatedRows]
# 
# #  save these lists as text in a special subfoler
# file1 <- here("Data","P30_P90-Grant","Shared_P30P90.txt")
# file2 <- here("Data","P30_P90-Grant","P30only.txt")
# file3 <- here("Data","P30_P90-Grant","P90only.txt")
# 
# write.table(x = Shared_P30P90, file1, sep = "\t")
# write.table(x = OnlyP30Grn, file2, sep = "\t")
# write.table(x = OnlyP90Grn, file3, sep = "\t")

# =======================================
# =======================================
# =======================================

