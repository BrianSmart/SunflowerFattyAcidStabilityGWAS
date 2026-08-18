# Genotype analysis, including PCA and clustering

library(pegas)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggrepel)
library(gridExtra)
library(dplyr)
library(AGHmatrix)
library(vcfR)
library(cluster)
library(factoextra)
library(dendextend)
library(ggfortify)
library(ggpubr)
library(poppr)
library(circlize)

setwd("~/Genotype_Analysis")

# create dosage matrix from VCF
# SAM085 & SAM272 missing
VCF <- read.vcfR("IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf")

# extract read depth per SNP. There is no depth column. Due to beagle imputation?
# dp <- extract.gt(VCF, element = "DP", as.numeric=TRUE)
# sum(is.na(dp))
# hist(na.omit(dp))

# extract and plot MAFs
maf <- maf(VCF)
hist(maf[,4], breaks=20, xlim = c(0,0.5))

# extract allele dosage from VCF
dosage <- extract.gt(VCF, as.numeric = F)

# replace VCF notation to numbers. Not elegant but works
dosage <- replace(dosage, dosage == "0|0", 0)
dosage <- replace(dosage, dosage == "0|1", 1)
dosage <- replace(dosage, dosage == "1|0", 1)
dosage <- replace(dosage, dosage == "1|1", 2)


dosage <- as.numeric(dosage)
# create correctly formatted matrix
matrix <- t(matrix(unlist(dosage), nrow = length(dosage)/287, ncol = 287))

# hierarchical clustering on genotype data
#nei_dist <- nei.dist(matrix) # calculate Nei's distance
euc_dist <- dist(matrix) # calculate Euclidean distance as alternative

#nei_dist_wardd2 <- (hclust(d= nei_dist, method = "ward.D2")) # cluster using ward.d2
euc_dist_wardd2 <- (hclust(d= euc_dist, method = "ward.D2")) # cluster using ward.d2


#nei_dend <- as.dendrogram(nei_dist_wardd2) # cluster to dendrogram
euc_dend <- as.dendrogram(euc_dist_wardd2) # cluster to dendrogram

#plot(nei_dend) # visualize dendrogram

# load file with FAD presence/absence data, color labels accordingly
fad <- read.csv("../Phenotype_Analysis/FAD2-1_mutant_status_renamed.csv")
labels_colors(euc_dend) <- fad$FAD2.1_flag[order.dendrogram(euc_dend)]

png(filename = "Euclidean_Dendrogram.png", width = 4500, height = 1000, units = "px", res = 200)
plot(euc_dend, main="Dendrogram using Euclidean distance") 
dev.off()

circlize_dendrogram(euc_dend,
                    labels_track_height = NA,
                    dend_track_height = 0.5)

# distance matrix visualization using euclidean distance.
png(filename = "distance_matrix.png", width = 2000, height = 2000, units = "px", res = 150)
fviz_dist(euc_dist) # doesn't look like population structure or duplicates
dev.off()

sum(is.infinite(matrix))
sum(is.na(matrix))

pca = prcomp(na.omit(matrix), scale. = T)
vars = 100*pca$sdev / sum(pca$sdev)

pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

png(filename = "pca_scree_plot.png", width = 1000, height = 1000, units = "px", res = 100)
barplot(pca.var.per, xlab="Principal Component", ylab="Percentage of Variation", cex.names = 3)
dev.off()
# 
# # extract PCs as covariates for GWAS
# # scrapped for now, vcf2gwas can do this automatically
# PCs <- get_pca(pca, element = "ind")
# covariates <- matrix(data=1, nrow=287, ncol=1)
# covariates <- cbind(covariates, PCs$coord[,1:3])
# 
# write.csv()

pc1 <- autoplot(pca, label = F,label.repel=F)+ theme_bw() #+ scale_color_manual(values=cbbPalette) #+ scale_color_brewer(palette="Set2")
pc2 <- autoplot(pca, x=2, y=3, label = F,label.repel=F)+ theme_bw() #+ scale_color_manual(values=cbbPalette)

png(filename = "PCA_with_scaling.png", width = 1250, height = 2000, units = "px", res = 200)
ggarrange(pc1,pc2,ncol = 1)
dev.off()
