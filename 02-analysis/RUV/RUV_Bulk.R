#Author Christine Muheim, code was adapet from the code in Risso et al., 2015 and Ingiosi et al., 2019.

library(limma)
library(edgeR)
library(EDASeq)
library(RUVSeq)
library(ffpe)
library(RColorBrewer)
library(DESeq2)
library(org.Mm.eg.db)
library(NOISeq)
library(stringr)
library(here)

sink('sessionInfo_RUV.txt')
sessionInfo()
sink()

# I am using projects to keep my data better organized
# ==============
# the code below uses the gene level summariedexperiemtn file created with txmeta in the previous step.
# cm, March 2022, I use the genome data first
# file <- here("Data","gse_mouse_sleep_complete.rds")
# gse_mouse_SleepDevelopment<- readRDS(file)
# gse_HCSD_WT_SleepDevelopment <- assays(gse_mouse_SleepDevelopment)[["counts"]]
# # ### save  as txt file
# file <- (here("Data","gse_HCSD_WT_SleepDevelopment_salmon.txt"))
# write.table(x = gse_HCSD_WT_SleepDevelopment, file = (file), sep = "\t")

# # note the gene-list contains version numbers after ".", remove everything with and after the "." on the gene level, this is not needed
# # open the tables

# file <- (here("Data","gse_HCSD_WT_SleepDevelopment_salmon.txt"))
# SDHC_wt_all <- read.table(file, row.names = 1, header = TRUE)
# colnames(SDHC_wt_all)
# rownames(SDHC_wt_all ) <- lapply(rownames(SDHC_wt_all),  sub, pattern = "\\.\\d+$", replacement = "")
# #  save again, indicated in the file name that this is now different!
# file <- (here("Data","gse_HCSD_WT_SleepDevelopment_salmon_trunc.txt"))
# write.table(x = SDHC_wt_all, file = (file), sep = "\t")

#  ==============
# once you have the text file made, load it into R
file <- (here("Data","gse_HCSD_WT_SleepDevelopment_salmon_trunc.txt"))
sdrs <- read.table(file, row.names = 1, header = TRUE)
# sdrs <- read.table("se_SD_HC_WT_SH3_P24_P30_ALLdata_salmon.txt", row.names = 1, header = TRUE)  #for the transcriptome data
head(sdrs)
colnames(sdrs)


# Change the column names for better visualization in the plots
# Here is how I did this
colnames(sdrs) <- c("HC16.1","HC16.2","HC16.3","HC16.4","HC16.5","HC24.1","HC24.2","HC24.3","HC24.4","HC24.5","HC24.6","HC24.7","HC24.8","HC30.1","HC30.2","HC30.3", "HC30.4","HC30.5","HC90.1","HC90.2","HC90.3","HC90.4","HC90.5","HC90.6","HC90.7","HC90.8","SD16.1","SD16.2","SD16.3","SD16.4","SD16.5","SD24.1","SD24.2","SD24.3","SD24.4","SD24.5","SD24.6","SD24.7","SD24.8","SD30.1","SD30.2","SD30.3","SD30.4","SD30.5","SD90.1","SD90.2","SD90.3","SD90.4","SD90.5","SD90.6","SD90.7","SD90.8")

colnames(sdrs)
dim(sdrs)
# [1] 54307    52

# Next, read the the positive control file. All of these files are published in the the eLife paper (Ingiosi et al.,2019)
sd.pos.ctrls <- read.table("PosControls_genes.txt")
dim(sd.pos.ctrls)
# [1] 572   1
x <- as.factor(rep(c("HC16","HC24","HC30","HC90","SD16","SD24","SD30","SD90"), c(5,8,5,8,5,8,5,8)))
names(x) <- colnames(sdrs)

# Filter out non expressed genes: We filter out any that are present less than 10x across 5 or more samples
# we choose to use a less stringent filtering approach to make sure do not bias towards genes expressed in the adult
filter <- apply(sdrs, 1, function(x) length(x[which(x > 10)]) > 5) 
filtered <- as.matrix(sdrs)[filter, ]
dim(filtered)
# [1] 18132    52

# ###### save the "all genes list", will be needed as background for DAVID
file <-here("Data","GeneListALLages_Filtered.txt")
write.table(x = filtered, file = file, sep = "\t")


# Use negative controls, as published in Ingiosi et al., 2019
sd.neg.ctrls <- read.table("NegControls_genes.txt")
dim(sd.neg.ctrls)
# [1] 6418    1

int.neg.ctrls <- sd.neg.ctrls
int.neg.ctrls <- unique(int.neg.ctrls)

colors <- brewer.pal(8, "Dark2")
colors2 <- brewer.pal(9,"Set1") 
colLib <- colors[x]

colLib2<- unique(colLib)
pch_vec<- c(15,16,17,18,19,20)

# UQ normalization. We tried TMM as well and found better recovery of positive controls with UQ.
uq <- betweenLaneNormalization(filtered, which = "upper")
dim(uq)

# Generating RLE Plots: RLE plots reveal confounders when the mean and the variance are not similar.
plot1 <- plotRLE(uq, col= colLib, outline = FALSE, las = 3, ylim = c(-0.7, 0.7), ylab = "Relative Log Expression",main="Upper Quantile", cex.axis = 1, cex.lab = 1)

plot2<- plotPCA(uq, labels=FALSE, pch=19, col = colLib,main="Upper Quantile", cex = 1.25, cex.axis = 1, cex.lab = 1, xlim = c(-0.3, 0.3), ylim = c(-0.43, 0.43))
# Add a legend
legend("topright", 
       legend = c("HC-P16","HC-P24","HC-P30","HC-P90","SD-P16","SD-P24","SD-P30","SD-P90"), 
       col = colLib2, 
       pch = 19, 
       bty = "n", 
       pt.cex = 1.4, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.025, 0.025))

# RUV Normalization with varying k. 
groups <- matrix(data = c(1:5,rep(-1,3),6:13,14:18,rep(-1,3),19:26,27:31,rep(-1,3),32:39,40:44,rep(-1,3),45:52), nrow = 8, byrow = TRUE)
# neg <- row.names(filtered)   # if you would rather use ALL of the genes rather than this specific list of negative controls
neg <- intersect(int.neg.ctrls[,1], rownames(uq))
s <- RUVs(x=uq, cIdx=neg, scIdx=groups, k=3)

# Following RUV Normalization, plot RLE and PCA: 
plot3 <- plotRLE(s$normalizedCounts, col = colLib, outline = FALSE, las = 3, ylim = c(-0.6, 0.6), ylab= "Relative Log Expression",main="RUV k=3", cex.axis = 1, cex.lab = 1)

plot4 <- plotPCA(s$normalizedCounts,labels = FALSE, pch = 19, col = colLib, cex = 1.25, cex.axis = 1, cex.lab = 1, main="RUV k=3",xlim = c(-0.3, 0.3), ylim = c(-0.3, 0.3))
# Add a legend
legend("topright",
       legend = c("HC-P16","HC-P24","HC-P30","HC-P90","SD-P16","SD-P24","SD-P30","SD-P90"),
       col = colLib2,
       pch = 19,
       bty = "n",
       pt.cex = 1.4,
       cex = 1.2,
       text.col = "black",
       horiz = F ,
       inset = c(0.025, 0.025))

# # positive controls(will be red in the volcano plots)
sd.pos <- intersect(sd.pos.ctrls[,1], rownames(filtered))
# length(sd.pos)
# # [1] 571 P30

 #   UPPER Quantile PLOTS, P16
# design <- model.matrix(~x)
# y <- DGEList(counts = filtered, group = x)
# y <- calcNormFactors(y, method = "upperquartile")
# y <- estimateDisp(y, design, verbose = TRUE)
# fit <- glmQLFit(y, design, robust = TRUE)
# qlf <- glmQLFTest(fit, coef=5) # compare HCP16 to SDP16, CORRECT way to select comparisons
# # #
#  summary(decideTests(qlf))
# topUQSD <- topTags(qlf, n = Inf)$table
# # ## Histogram for UQ
# hist(topUQSD$PValue, main = "", xlab = "p-value", breaks = 100, ylim = c(0, 2250))
# ## Volcano Plot
# plot(topUQSD[, 1], -log10(topUQSD$PValue), pch = 20, col = "gray", cex= 0.5, ylab = "-log10(p-value)", xlab= "log2(SD/HC)",main="Upper Quantile P16", ylim = c(0, 15), xlim = c(-3, 3), cex.lab = 1, cex.axis = 1)
# # de is blue
# de <- rownames(topUQSD[topUQSD$FDR <= 0.05, ]) # In eLife it was less than 0.05, although I noticed some of Lucia's old papers used 0.01. I stuck with 0.05 to keep it similar to eLife
# points(topUQSD[de, 1], -log10(topUQSD[de, "PValue"]), pch = 20, col = colors2[2], cex = 0.5, lwd = 2)
#
# points(topUQSD[sd.pos, 1], -log10(topUQSD[sd.pos, "PValue"]), pch = 20, col = colors2[1], cex = 0.5, lwd = 2)
#
# (sum(topUQSD[sd.pos, "FDR"] < .05))/(length(sd.pos) )* 100
# (sum(topUQSD[sd.pos, "FDR"] < .05))
# difex1 <- sum(topUQSD$FDR < 0.05)

# #   UPPER Quantile PLOTS, P24
# design <- model.matrix(~x)
# y <- DGEList(counts = filtered, group = x)
# y <- calcNormFactors(y, method = "upperquartile")
# y <- estimateDisp(y, design, verbose = TRUE)
# fit <- glmQLFit(y, design, robust = TRUE)
# qlf <- glmQLFTest(fit, contrast=c(0,-1,0,0,0,1,0,0)) # compare HCP24 to SDP24,
# #
# summary(decideTests(qlf))
# topUQSD <- topTags(qlf, n = Inf)$table
# # ## Histogram for UQ
# hist(topUQSD$PValue, main = "", xlab = "p-value", breaks = 100, ylim = c(0, 2250))
# ## Volcano Plot
# plot(topUQSD[, 1], -log10(topUQSD$PValue), pch = 20, col = "gray", cex= 0.5, ylab = "-log10(p-value)", xlab= "log2(SD/HC)", main="Upper Quantile P24", ylim = c(0, 15), xlim = c(-3, 3), cex.lab = 1, cex.axis = 1)
# # de is blue
# de <- rownames(topUQSD[topUQSD$FDR <= 0.05, ]) # In eLife it was less than 0.05, although I noticed some of Lucia's old papers used 0.01. I stuck with 0.05 to keep it similar to eLife
# points(topUQSD[de, 1], -log10(topUQSD[de, "PValue"]), pch = 20, col = colors2[2], cex = 0.5, lwd = 2)
#
# points(topUQSD[sd.pos, 1], -log10(topUQSD[sd.pos, "PValue"]), pch = 20, col = colors2[1], cex = 0.5, lwd = 2)
#
# # negative controls are green. If you want to plot them you can, however Lucia and Davide think it is redundant and not necessary
# # points(topUQSD[neg, 1], -log10(topUQSD[neg, "PValue"]), pch = 1, col = colors[3], cex = 1, lwd = 2)
#
# (sum(topUQSD[sd.pos, "FDR"] < .05))/(length(sd.pos) )* 100
# (sum(topUQSD[sd.pos, "FDR"] < .05))
# difex1 <- sum(topUQSD$FDR < 0.05)
#alternative use length(de)
# #
#
# #   UPPER Quantile PLOTS, P30
# design <- model.matrix(~x)
# y <- DGEList(counts = filtered, group = x)
# y <- calcNormFactors(y, method = "upperquartile")
# y <- estimateDisp(y, design, verbose = TRUE)
# fit <- glmQLFit(y, design, robust = TRUE)
# qlf <- glmQLFTest(fit, contrast=c(0,0,-1,0,0,0,1,0))   # compare HCP30 to SDP30
# 
# summary(decideTests(qlf))
# topUQSD <- topTags(qlf, n = Inf)$table
# ## Histogram for UQ
# hist(topUQSD$PValue, main = "", xlab = "p-value", breaks = 100, ylim = c(0, 2250))
# ## Volcano Plot
# plot(topUQSD[, 1], -log10(topUQSD$PValue), pch = 20, col = "gray", cex= 0.5, ylab = "-log10(p-value)", xlab= "log2(SD/HC)",main="Upper Quantile P30", ylim = c(0, 15), xlim = c(-3, 3), cex.lab = 1, cex.axis = 1)
# # de is blue
# de <- rownames(topUQSD[topUQSD$FDR <= 0.05, ]) # In eLife it was less than 0.05, although I noticed some of Lucia's old papers used 0.01. I stuck with 0.05 to keep it similar to eLife
# points(topUQSD[de, 1], -log10(topUQSD[de, "PValue"]), pch = 20, col = colors2[2], cex = 0.5, lwd = 2)
#
# points(topUQSD[sd.pos, 1], -log10(topUQSD[sd.pos, "PValue"]), pch = 20, col = colors2[1], cex = 0.5, lwd = 2)
#
# (sum(topUQSD[sd.pos, "FDR"] < .05))/(length(sd.pos) )* 100
# (sum(topUQSD[sd.pos, "FDR"] < .05))
# difex1 <- sum(topUQSD$FDR < 0.05)
# #alternative use length(de)

# #   UPPER Quantile PLOTS, P90
# design <- model.matrix(~x)
# y <- DGEList(counts = filtered, group = x)
# y <- calcNormFactors(y, method = "upperquartile")
# y <- estimateDisp(y, design, verbose = TRUE)
# fit <- glmQLFit(y, design, robust = TRUE)
# qlf <- glmQLFTest(fit, contrast=c(0,0,0,-1,0,0,0,1))   # compare HCP90 to SDP90#
# 
# summary(decideTests(qlf))
# topUQSD <- topTags(qlf, n = Inf)$table
# ## Histogram for UQ
# hist(topUQSD$PValue, main = "", xlab = "p-value", breaks = 100, ylim = c(0, 2250))
# ## Volcano Plot
# plot(topUQSD[, 1], -log10(topUQSD$PValue), pch = 20, col = "gray", cex= 0.5, ylab = "-log10(p-value)", xlab= "log2(SD/HC)",main="Upper Quantile P90", ylim = c(0, 15), xlim = c(-3, 3), cex.lab = 1, cex.axis = 1)
# # de is blue
# de <- rownames(topUQSD[topUQSD$FDR <= 0.05, ]) # In eLife it was less than 0.05, although I noticed some of Lucia's old papers used 0.01. I stuck with 0.05 to keep it similar to eLife
# points(topUQSD[de, 1], -log10(topUQSD[de, "PValue"]), pch = 20, col = colors2[2], cex = 0.5, lwd = 2)
#
# points(topUQSD[sd.pos, 1], -log10(topUQSD[sd.pos, "PValue"]), pch = 20, col = colors2[1], cex = 0.5, lwd = 2)
#
# (sum(topUQSD[sd.pos, "FDR"] < .05))/(length(sd.pos) )* 100
# (sum(topUQSD[sd.pos, "FDR"] < .05))
# difex1 <- sum(topUQSD$FDR < 0.05)
# #alternative use length(de)

# ==========================================
### Differential Expression: RUV Normalization P16

design <- model.matrix(~x + s$W)
y1 <- DGEList(counts = filtered, group = x)
y1 <- calcNormFactors(y1, method = "upperquartile")
y1 <- estimateDisp(y1, design, verbose = TRUE)
fit1 <- glmQLFit(y1, design, robust = TRUE)
qlf1 <- glmQLFTest(fit1, coef=5) # compare HCP16 to SDP16, this one has to done like this, it compares to the first group

summary(decideTests(qlf1))
topSD16 <- topTags(qlf1, n = Inf)$table

sd.pos <- intersect(sd.pos.ctrls[,1], rownames(filtered))
length(sd.pos)
## Determine the percentage of positive control genes detected
(sum(topSD16[sd.pos, "FDR"] < .05))/(length(sd.pos)) * 100
(sum(topSD16[sd.pos, "FDR"] < .05))

# #  export all the DE genes for each comparison as txt file
de <- rownames(topSD16[topSD16$FDR < 0.05, ]) # keep this at 0.05 for the significant genes, also export all genes further below
DE16genes <- (topSD16[de, ] )
dim(DE16genes)
# 
# # save all the gens with FDR, logFC, logCPM, F and Pvalue
file1 <-here("Data","HCvsSD_P16.txt")
write.table(x = DE16genes, file = file1, sep = "\t")
# # save the entire list, just in case
#deALL <- rownames(topSD16[topSD16$FDR <= 1, ]) # might not be necessary to run again, changed this to 1, if you want the signifcants only change it back to 0.05 or less
#DE16genesALL <- (topSD16[deALL, ] )
#file2 <-here("Data","HCvsSD_P16_ALL.txt")
#write.table(x = DE16genesALL, file = file2, sep = "\t")

# Histogram
plot5 <- hist(topSD16$PValue, main = "P16", xlab= "p-value", breaks = 100, ylim = c(0, 2200))

## Volcano Plot
plot6 <- plot(topSD16[, 1], -log10(topSD16$PValue), pch = 20, col = "gray", cex = 0.5, ylab = "-log10(p-value)", xlab = "log2(SD/HC)",main="P16", ylim = c(0, 25), xlim = c(-4, 4), cex.lab = 1, cex.axis = 1)

# de genes are blue
de <- rownames(topSD16[topSD16$FDR <= 0.05, ])
points(topSD16[de, 1], -log10(topSD16[de, "PValue"]), pch = 20, col = colors2[2], cex = 0.5)
points(topSD16[de, 1], -log10(topSD16[de, "PValue"]), pch = 20, col = colors2[2], cex = 0.5)

# positive controls are red
points(topSD16[sd.pos, 1], -log10(topSD16[sd.pos, "PValue"]), pch = 20, col = colors2[1], cex =0.5, lwd = 2)

# negative controls are green
# points(topSD[neg, 1], -log10(topSD[neg, "PValue"]), pch = 1, col = colors[3], cex = 1, lwd = 2)


### Differential Expression: RUV Normalization P24
design <- model.matrix(~x + s$W)
y1 <- DGEList(counts = filtered, group = x)
y1 <- calcNormFactors(y1, method = "upperquartile")
y1 <- estimateDisp(y1, design, verbose = TRUE)
fit2 <- glmQLFit(y1, design, robust = TRUE)
qlf2 <- glmQLFTest(fit2, contrast=c(0,-1,0,0,0,1,0,0,0,0,0)) # compare HCP24 to SDP24, 

summary(decideTests(qlf2))
topSD24 <- topTags(qlf2, n = Inf)$table

sd.pos <- intersect(sd.pos.ctrls[,1], rownames(filtered))
length(sd.pos)
## Determine the percentage of positive control genes detected
(sum(topSD24[sd.pos, "FDR"] < .05))/(length(sd.pos)) * 100
(sum(topSD24[sd.pos, "FDR"] < .05))

# #  export all the DE genes for each comparison as txt file
de <- rownames(topSD24[topSD24$FDR < 0.05, ]) # keep this at 0.05 for the significant genes, also export all genes further below
DE24genes <- (topSD24[de, ] )
dim(DE24genes)
# #  save all the gens with FDR, logFC, logCPM, F and Pvalue!
file1 <-here("Data","HCvsSD_P24.txt")
write.table(x = DE24genes, file = file1, sep = "\t")
# 
# # # save the entire list too, just in case
#deALL <- rownames(topSD24[topSD24$FDR <= 1, ]) # might not be necessary to run again, changed this to 1, if you want the signifcants only change it back to 0.05 or less
#DE24genesALL <- (topSD24[deALL, ] )
#file2 <-here("Data","HCvsSD_P24_ALL.txt")
#write.table(x = DE24genesALL, file = file2, sep = "\t")

# Histogram
plot7 <- hist(topSD24$PValue, main = "P24", xlab= "p-value", breaks = 100, ylim = c(0, 2200))

## Volcano Plot
plot8 <- plot(topSD24[, 1], -log10(topSD24$PValue), pch = 20, col = "gray", cex = 0.5, ylab = "-log10(p-value)", xlab = "log2(SD/HC)",main="P24", ylim = c(0, 25), xlim = c(-5, 5), cex.lab = 1, cex.axis = 1)

# de genes are blue
de <- rownames(topSD24[topSD24$FDR <= 0.05, ])
points(topSD24[de, 1], -log10(topSD24[de, "PValue"]), pch = 20, col = colors2[2], cex = 0.5)
points(topSD24[de, 1], -log10(topSD24[de, "PValue"]), pch = 20, col = colors2[2], cex = 0.5)

# positive controls are red
points(topSD24[sd.pos, 1], -log10(topSD24[sd.pos, "PValue"]), pch = 20, col = colors2[1], cex =0.5, lwd = 2)

# negative controls are green
# points(topSD[neg, 1], -log10(topSD[neg, "PValue"]), pch = 1, col = colors[3], cex = 1, lwd = 2)




### Differential Expression: RUV Normalization P30
design <- model.matrix(~x + s$W)
y1 <- DGEList(counts = filtered, group = x)
y1 <- calcNormFactors(y1, method = "upperquartile")
y1 <- estimateDisp(y1, design, verbose = TRUE)
fit3 <- glmQLFit(y1, design, robust = TRUE)
qlf3 <- glmQLFTest(fit3, contrast=c(0,0,-1,0,0,0,1,0,0,0,0))   # compare HCP30 to SDP30

summary(decideTests(qlf3))
topSD30 <- topTags(qlf3, n = Inf)$table

sd.pos <- intersect(sd.pos.ctrls[,1], rownames(filtered))
length(sd.pos)
## Determine the percentage of positive control genes detected
(sum(topSD30[sd.pos, "FDR"] < .05))/(length(sd.pos)) * 100
(sum(topSD30[sd.pos, "FDR"] < .05))

# Histogram
plot9 <- hist(topSD30$PValue, main = "P30", xlab= "p-value", breaks = 100, ylim = c(0, 2200))

## Volcano Plot
plot10 <- plot(topSD30[, 1], -log10(topSD30$PValue), pch = 20, col = "gray", cex = 0.5, ylab = "-log10(p-value)", xlab = "log2(SD/HC)",main="P30", ylim = c(0, 25), xlim = c(-6, 6), cex.lab = 1, cex.axis = 1)

# de genes are blue
de <- rownames(topSD30[topSD30$FDR <= 0.05, ])
points(topSD30[de, 1], -log10(topSD30[de, "PValue"]), pch = 20, col = colors2[2], cex = 0.5)
points(topSD30[de, 1], -log10(topSD30[de, "PValue"]), pch = 20, col = colors2[2], cex = 0.5)

# positive controls are red
points(topSD30[sd.pos, 1], -log10(topSD30[sd.pos, "PValue"]), pch = 20, col = colors2[1], cex =0.5, lwd = 2)

# negative controls are green
# points(topSD[neg, 1], -log10(topSD[neg, "PValue"]), pch = 1, col = colors[3], cex = 1, lwd = 2)

# #  export all the DE genes for each comparison as txt file
de <- rownames(topSD30[topSD30$FDR < 0.05, ]) # keep this at 0.05 for the significants, also export all genes further below
DE30genes <- (topSD30[de, ] )
dim(DE30genes)
# # save all the gens with FDR, logFC, logCPM, F and Pvalue!
file1 <-here("Data","HCvsSD_P30.txt")
write.table(x = DE30genes, file = file1, sep = "\t")

# # save the entire list too, just in case
#deALL <- rownames(topSD30[topSD30$FDR <= 1, ]) # might not be necessary to run again, changed this to 1, if you want the signifcants only change it back to 0.05 or less
#DE30genesALL <- (topSD30[deALL, ] )
#file2 <-here("Data","HCvsSD_P30_ALL.txt")
#write.table(x = DE30genesALL, file = file2, sep = "\t")

### Differential Expression: RUV Normalization P90
design <- model.matrix(~x + s$W)
y1 <- DGEList(counts = filtered, group = x)
y1 <- calcNormFactors(y1, method = "upperquartile")
y1 <- estimateDisp(y1, design, verbose = TRUE)
fit4 <- glmQLFit(y1, design, robust = TRUE)
qlf4 <- glmQLFTest(fit4, contrast=c(0,0,0,-1,0,0,0,1,0,0,0))   # compare HCP90 to SDP90#
# qlf4 <- glmQLFTest(fit4, coef=4) # compare HCP90 to SDP90, alternative way to selectect comparisons
summary(decideTests(qlf4))

topSD90 <- topTags(qlf4, n = Inf)$table
sd.pos <- intersect(sd.pos.ctrls[,1], rownames(filtered))
length(sd.pos)
## Determine the percentage of positive control genes detected
(sum(topSD90[sd.pos, "FDR"] < .05))/(length(sd.pos)) * 100
(sum(topSD90[sd.pos, "FDR"] < .05))

# Histogram
plot11 <- hist(topSD90$PValue, main = "P90", xlab= "p-value", breaks = 100, ylim = c(0, 2200))

## Volcano Plot
plot12 <- plot(topSD90[, 1], -log10(topSD90$PValue), pch = 20, col = "gray", cex = 0.5, ylab = "-log10(p-value)", xlab = "log2(SD/HC)",main="P90", ylim = c(0, 25), xlim = c(-7, 7), cex.lab = 1, cex.axis = 1)

# de genes are blue
de <- rownames(topSD90[topSD90$FDR <= 0.05, ])
points(topSD90[de, 1], -log10(topSD90[de, "PValue"]), pch = 20, col = colors2[2], cex = 0.5)
points(topSD90[de, 1], -log10(topSD90[de, "PValue"]), pch = 20, col = colors2[2], cex = 0.5)

# positive controls are red
points(topSD90[sd.pos, 1], -log10(topSD90[sd.pos, "PValue"]), pch = 20, col = colors2[1], cex =0.5, lwd = 2)

# negative controls are green
# points(topSD[neg, 1], -log10(topSD[neg, "PValue"]), pch = 1, col = colors[3], cex = 1, lwd = 2)

# #  export all the DE genes for each comparison as txt file
de <- rownames(topSD90[topSD90$FDR < 0.05, ]) # keep this at 0.05 for the significants, also export all genes further below
DE90genes <- (topSD90[de, ] )
dim(DE90genes)
# #  save all the gens with FDR, logFC, logCPM, F and Pvalue!
file1 <-here("Data","HCvsSD_P90.txt")
write.table(x = DE90genes, file = file1, sep = "\t")
# 
# # save the entire list too, just in case
#deALL <- rownames(topSD90[topSD90$FDR <= 1, ]) # might not be necessary to run again, changed this to 1, if you want the signifcants only change it back to 0.05 or less
#DE90genesALL <- (topSD90[deALL, ] )
#file2 <-here("Data","HCvsSD_P90_ALL.txt")
#write.table(x = DE90genesALL, file = file2, sep = "\t")

