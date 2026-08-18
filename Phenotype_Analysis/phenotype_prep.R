# phenotype analysis and formatting for GWAS

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggfortify)
library(cluster)
library(factoextra)
library(dendextend)
library(readxl)
library(dplyr)
library(pheatmap) 
library(vegan)
library(stringr)
library(cellWise)
library(reshape)
library(gge)
library(PerformanceAnalytics)
library(hrbrthemes)
library(GGally)
library(viridis)

setwd("D:/Sunflower_stability_GWAS/Phenotype_Analysis/")

# colorblind color palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# import phenotypes from separate file, as the source data was updated

MN_2015 <- as.data.frame(read_excel("MN15 mean.xlsx", sheet = 1, col_types = c("text", rep("numeric", 8)))[c(-1,-272),]) 
MN_2015$PPN <- paste0("SAM", str_pad(as.character(MN_2015$PPN),3,pad = "0"))
rownames(MN_2015) <- MN_2015$PPN
colnames(MN_2015) <- c("Line", "Palm", "Stea", "Oleic", "Lino","n1", "n2", "n3", "n4")
#MN_2015[c(4,69,116,255),2:9]


MN_2016_day_1 <- as.data.frame(read_excel("MN2016early mean stderr.xlsx", sheet = 1, col_types = c("text", rep("numeric", 8)))[c(-1,-272),]) 
MN_2016_day_1$PPN <- paste0("SAM", str_pad(as.character(MN_2016_day_1$PPN),3,pad = "0"))
rownames(MN_2016_day_1) <- MN_2016_day_1$PPN
colnames(MN_2016_day_1) <- c("Line", "Palm", "Stea", "Oleic", "Lino","n1", "n2", "n3", "n4")
#MN_2016_day_1[c(4,69,116,255),2:9] # check if questionable genotypes were removed prior

MN_2016_day_2 <- as.data.frame(read_excel("MN2016late mean stderr.xlsx", sheet = 1, col_types = c("text", rep("numeric", 8)))[c(-1,-272),]) 
MN_2016_day_2$PPN <- paste0("SAM", str_pad(as.character(MN_2016_day_2$PPN),3,pad = "0"))
rownames(MN_2016_day_2) <- MN_2016_day_2$PPN
colnames(MN_2016_day_2) <- c("Line", "Palm", "Stea", "Oleic", "Lino","n1", "n2", "n3", "n4")
#MN_2016_day_2[c(4,69,116,255),2:9]

GA <- as.data.frame(read_excel("GA mean stderr.xlsx", sheet = 1, col_types = c("text", rep("numeric", 8)))[c(-1,-272),]) 
GA$PPN <- paste0("SAM", str_pad(as.character(GA$PPN),3,pad = "0"))
rownames(GA) <- GA$PPN
colnames(GA) <- c("Line", "Palm", "Stea", "Oleic", "Lino","n1", "n2", "n3", "n4")
#GA[c(4,69,116,255),2:9]

BC <- as.data.frame(read_excel("BC mean stderr.xlsx", sheet = 1, col_types = c("text", rep("numeric", 8)))[c(-1,-272),]) 
BC$PPN <- paste0("SAM", str_pad(as.character(BC$PPN),3,pad = "0"))
rownames(BC) <- BC$PPN
colnames(BC) <- c("Line", "Palm", "Stea", "Oleic", "Lino","n1", "n2", "n3", "n4")
#BC[c(4,69,116,255),2:9]


IA <- as.data.frame(read_excel("IA mean.xlsx", sheet = 1, col_types = c("text", rep("numeric", 8)))[c(-1,-272),]) 
IA$PPN <- paste0("SAM", str_pad(as.character(IA$PPN),3,pad = "0"))
rownames(IA) <- IA$PPN
colnames(IA) <- c("Line", "Palm", "Stea", "Oleic", "Lino","n1", "n2", "n3", "n4")
IA[c(4,69,116,255),2:9]

IA_2013 <- as.data.frame(read_excel("IA 2013-2014 means stderr.xlsx", sheet = 1, col_types = c("text", rep("numeric", 8)))[c(-1,-272),]) # 188, 252 were missing, added them manually in excel
IA_2013$PPN <- paste0("SAM", str_pad(as.character(IA_2013$PPN),3,pad = "0")) # add "SAM" in front of each individual's number
rownames(IA_2013) <- IA_2013$PPN
colnames(IA_2013) <- c("Line", "Palm", "Stea", "Oleic", "Lino","n1", "n2", "n3", "n4")
IA_2013[c(4,69,116),2:9] <- NA # remove questionable genotypes

IA_2014 <- as.data.frame(read_excel("IA 2013-2014 means stderr.xlsx", sheet = 2, col_types = c("text", rep("numeric", 8)))[c(-1,-272),]) # 188, 252 were missing, added them manually in excel
IA_2014$PPN <- paste0("SAM", str_pad(as.character(IA_2014$PPN),3,pad = "0")) # add "SAM" in front of each individual's number
rownames(IA_2014) <- IA_2014$PPN
colnames(IA_2014) <- c("Line", "Palm", "Stea", "Oleic", "Lino","n1", "n2", "n3", "n4")
IA_2014[c(4,69,116),2:9] <- NA # remove questionable genotypes


# replace phenotypes based on only one sample with NA, as these might be unreliable

MN_2015[!is.na(MN_2015$n1) & MN_2015$n1=="1", c(2:5)] <- NA 
MN_2016_day_1[!is.na(MN_2016_day_1$n1) & MN_2016_day_1$n1=="1", c(2:5)] <- NA 
MN_2016_day_2[!is.na(MN_2016_day_2$n1) & MN_2016_day_2$n1=="1", c(2:5)] <- NA 
GA[!is.na(GA$n1) & GA$n1=="1", c(2:5)] <- NA 
BC[!is.na(BC$n1) & BC$n1=="1", c(2:5)] <- NA 
IA[!is.na(IA$n1) & IA$n1=="1", c(2:5)] <- NA # only half left in this case
IA_2013[!is.na(IA_2013$n1) & IA_2013$n1=="1", c(2:5)] <- NA
IA_2014[!is.na(IA_2014$n1) & IA_2014$n1=="1", c(2:5)] <- NA

# remove n columns
MN_2015 <- MN_2015 [,c(1:5)]
MN_2016_day_1 <- MN_2016_day_1[,c(1:5)]
MN_2016_day_2 <- MN_2016_day_2[,c(1:5)]
BC <- BC[,c(1:5)]
GA <- GA[,c(1:5)]
IA <- IA[,c(1:5)]
IA_2013 <- IA_2013[,c(1:5)]
IA_2014 <- IA_2014[,c(1:5)]

# number of phenotypes left
sum(!is.na(MN_2015$Palm))
sum(!is.na(MN_2016_day_1$Palm))
sum(!is.na(MN_2016_day_2$Palm))
sum(!is.na(BC$Palm)) 
sum(!is.na(GA$Palm)) 
sum(!is.na(IA$Palm))
sum(!is.na(IA_2014$Palm))
sum(!is.na(IA_2013$Palm))

plot(IA_2014)
plot(IA)
#######################################################################
### comment out before writing to csv ###

# GGE biplot

BC$Region <- "BC" # add sample region as a column to easier differentiate them after merging
GA$Region <- "GA"
IA$Region <- "IA"
MN_2015$Region <- "MN_2015"
MN_2016_day_1$Region <- "MN_2016_day_1"
MN_2016_day_2$Region <- "MN_2016_day_2"
IA_2013$Region <- "IA_2013"
IA_2014$Region <- "IA_2014"

all_phenotypes <- bind_rows(BC,IA,GA,MN_2015,MN_2016_day_1,MN_2016_day_2, IA_2013, IA_2014)
#all_phenotypes <- bind_rows(BC,GA,MN_2016_day_1,MN_2016_day_2, IA_2013, IA_2014)

png("fatty_acid_relationships.png", width = 1000, height = 1000, res = 100)
#plot(all_phenotypes[-c(1,6)])
chart.Correlation(all_phenotypes[-c(1,6)])
dev.off()

cor.test(all_phenotypes$Oleic, all_phenotypes$Lino)
cor.test(all_phenotypes$Oleic, all_phenotypes$Stea)
cor.test(all_phenotypes$Oleic, all_phenotypes$Palm)


all_phenotypes_melt <- melt(all_phenotypes)

all_phenotypes_melt$ID <- with(all_phenotypes_melt, paste(Region, variable, sep = "_"))

gge <- gge(na.omit(all_phenotypes_melt), value ~ Line*ID, scale=T)

png("GGE_biplot.png", width = 1000, height = 1000, res = 120)
biplot(gge, origin = 0, hull=F)
dev.off()

png("GGE_biplot_addition.png", width = 1000, height = 1000, res = 120)
plot(gge)
dev.off()
# try different package
#library(metan)
# to-do
#####

BC$Line <- paste0((BC$Line),"_BC") # add location to each sample name
IA$Line <- paste0((IA$Line),"_IA")
GA$Line <- paste0((GA$Line),"_GA")
MN_2015$Line <- paste0((MN_2015$Line),"_MN_2015")
MN_2016_day_1$Line <- paste0((MN_2016_day_1$Line),"_MN_2016_day_1")
MN_2016_day_2$Line <- paste0((MN_2016_day_2$Line),"_MN_2016_day_2")
IA_2013$Line <- paste0((IA_2013$Line),"_IA_2013") # add location to each sample name
IA_2014$Line <- paste0((IA_2014$Line),"_IA_2014") # add location to each sample name


BC$Region <- "BC" # add sample region as a column to easier differentiate them after merging
GA$Region <- "GA"
IA$Region <- "IA"
MN_2015$Region <- "MN_2015"
MN_2016_day_1$Region <- "MN_2016_day_1"
MN_2016_day_2$Region <- "MN_2016_day_2"
IA_2013$Region <- "IA_2013"
IA_2014$Region <- "IA_2014"




# stack data condition-wise so that every oil is in one column
# doesn't work with identical rownames
phenotypes_for_GWAS <- bind_rows(BC,IA,GA,MN_2015,MN_2016_day_1,MN_2016_day_2, IA_2013, IA_2014)
rownames(phenotypes_for_GWAS) <- phenotypes_for_GWAS$Line

# PCA on fatty acid composition with all samples in all conditions
pca <- prcomp(na.omit(phenotypes_for_GWAS[,c(-1,-6)]), scale. = F)

vars = 100*pca$sdev / sum(pca$sdev)

pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
barplot(pca.var.per, xlab="Principal Component", ylab="Percentage of Variation", cex.names = 3)

png("PCA_phenotypes.png", width = 1000, height = 1000, res = 120)
autoplot(pca, label = T,label.repel=T, data = na.omit(phenotypes_for_GWAS), colour = 'Region') + 
  theme_bw() + scale_color_manual(values=cbbPalette) # outliers seem to coincide with low linoleic acid
dev.off()
#autoplot(pca, x=2, y=3, label = F,label.repel=F, data = na.omit(phenotypes_for_GWAS), colour = 'Region') + 
#  theme_bw() + scale_color_manual(values=cbbPalette)
# 
# # try hierarchical clustering on phenotypes
# dist <- dist(na.omit(phenotypes_for_GWAS))
# dist_wardd2 <- hclust(d=dist, method = "ward.D2")
# dend <- as.dendrogram(dist_wardd2)
# 
# dend %>% set("branches_k_color", value = cbbPalette, k = 6) %>% 
#   plot

# Pheatmaps to visualize phenotypes. 
# Trends of oil contents between samples persists for the most part, two exeptions are stearic in BC

png("BC_pheatmap.png", width = 1000, height = 1000, res = 120)
pheatmap(na.omit(BC[,-c(1,6)]), color = rev(hcl.colors(50, "Sunset")), cex = 1.5, cluster_cols = F, fontsize = 5, show_rownames = F)
dev.off()

pheatmap(na.omit(GA[,-c(1,6)]), color = rev(hcl.colors(50, "Sunset")), cex = 1.5, cluster_cols = F, fontsize = 5)
pheatmap(na.omit(IA[,-c(1,6)]), color = rev(hcl.colors(50, "Sunset")), cex = 1.5, cluster_cols = F, fontsize = 5)
pheatmap(na.omit(MN_2015[,-c(1,6)]), color = rev(hcl.colors(50, "Sunset")), cex = 1.5, cluster_cols = F, fontsize = 5)
pheatmap(na.omit(MN_2016_day_1[,-c(1,6)]), color = rev(hcl.colors(50, "Sunset")), cex = 1.5, cluster_cols = F, fontsize = 5)
pheatmap(na.omit(MN_2016_day_2[,-c(1,6)]), color = rev(hcl.colors(50, "Sunset")), cex = 1.5, cluster_cols = F, fontsize = 5)
pheatmap(na.omit(IA_2013[,-c(1,6)]), color = rev(hcl.colors(50, "Sunset")), cex = 1.5, cluster_cols = F, fontsize = 5)
pheatmap(na.omit(IA_2014[,-c(1,6)]), color = rev(hcl.colors(50, "Sunset")), cex = 1.5, cluster_cols = F, fontsize = 5)


pheatmap(na.omit(phenotypes_for_GWAS[,-c(1,6)]), color = rev(hcl.colors(50, "Sunset")), cex = 1, cluster_cols = F, fontsize = 5)
pheatmap(na.omit(phenotypes_for_GWAS[sample(length(phenotypes_for_GWAS$Palm), 200),-c(1,6)]), color = rev(hcl.colors(50, "Sunset")), cex = 1, cluster_cols = F, fontsize = 5)
#######################################################################

# write files to csv
write.table(BC, file = "BC_oil_contents.csv", sep = ",", row.names = F, col.names = T)
write.table(IA, file = "IA_oil_contents.csv", sep = ",", row.names = F, col.names = T)
write.table(GA, file = "GA_oil_contents.csv", sep = ",", row.names = F, col.names = T)
write.table(MN_2015, file = "MN_2015_oil_contents.csv", sep = ",", row.names = F, col.names = T)
write.table(MN_2016_day_1, file = "MN_2016_day_1_oil_contents.csv", sep = ",", row.names = F, col.names = T)
write.table(MN_2016_day_2, file = "MN_2016_day_2_oil_contents.csv", sep = ",", row.names = F, col.names = T)
write.table(IA_2013, file = 'IA_2013_oil_contents.csv', sep = ",", row.names = F, col.names = T)
write.table(IA_2014, file = 'IA_2014_oil_contents.csv', sep = ",", row.names = F, col.names = T)
#write.table(phenotypes_for_GWAS, file = "All_locations_oil_contents.csv", sep = ",", row.names = F, col.names = T)

# (exclude FAD mutants). Set values to NA instead so they are ignored by vcf2GWAS
FAD_status <- read.csv("FAD2-1_mutant_status_renamed.csv")
#MN_2016_day_1_no_FAD <- MN_2016_day_1[FAD_status$FAD2.1_flag == 0,]
#MN_2016_day_1[FAD_status$FAD2.1_flag == 1,2:5] <- NA
IA_2014[FAD_status$FAD2.1_flag == 1,2:5] <- NA

# make a simple list of FAD2 mutants
fad <- FAD_status[FAD_status$FAD2.1_flag == 1,]$Line
write.table(fad, file = "fad_mutants_list.csv", sep = ",", row.names = F, col.names = F)


#write.table(MN_2016_day_1, file = "MN_2016_day_1_oil_contents_no_FAD.csv", sep = ",", row.names = F, col.names = T)
write.table(IA_2014, file = "IA_2014_oil_contents_no_FAD.csv", sep = ",", row.names = F, col.names = T)


################################################
### process standard errors in each location ###
################################################

MN_2016_day_1_stderr <- as.data.frame(read_excel("MN2016early mean stderr.xlsx", sheet = 2, col_types = c("text", rep("numeric", 8)))[c(-1,-272),c(1:5)])
MN_2016_day_1_stderr$PPN <- paste0("SAM", str_pad(as.character(MN_2016_day_1_stderr$PPN),3,pad = "0")) # add "SAM" in front of each individual's number
rownames(MN_2016_day_1_stderr) <- MN_2016_day_1_stderr$PPN
colnames(MN_2016_day_1_stderr) <- c("Line", "Palm_stderr", "Stea_stderr", "Oleic_stderr", "Lino_stderr")


MN_2016_day_2_stderr <- as.data.frame(read_excel("MN2016late mean stderr.xlsx", sheet = 2, col_types = c("text", rep("numeric", 8)))[c(-1,-272),c(1:5)])
MN_2016_day_2_stderr$PPN <- paste0("SAM", str_pad(as.character(MN_2016_day_2_stderr$PPN),3,pad = "0")) # add "SAM" in front of each individual's number
rownames(MN_2016_day_2_stderr) <- MN_2016_day_2_stderr$PPN
colnames(MN_2016_day_2_stderr) <- c("Line", "Palm_stderr", "Stea_stderr", "Oleic_stderr", "Lino_stderr")


BC_stderr <- as.data.frame(read_excel("BC mean stderr.xlsx", sheet = 2, col_types = c("text", rep("numeric", 8)))[c(-1,-272),c(1:5)])
BC_stderr$PPN <- paste0("SAM", str_pad(as.character(BC_stderr$PPN),3,pad = "0")) # add "SAM" in front of each individual's number
rownames(BC_stderr) <- BC_stderr$PPN
colnames(BC_stderr) <- c("Line", "Palm_stderr", "Stea_stderr", "Oleic_stderr", "Lino_stderr")


GA_stderr <- as.data.frame(read_excel("GA mean stderr.xlsx", sheet = 2, col_types = c("text", rep("numeric", 8)))[c(-1,-272),c(1:5)])
GA_stderr$PPN <- paste0("SAM", str_pad(as.character(GA_stderr$PPN),3,pad = "0")) # add "SAM" in front of each individual's number
rownames(GA_stderr) <- GA_stderr$PPN
colnames(GA_stderr) <- c("Line", "Palm_stderr", "Stea_stderr", "Oleic_stderr", "Lino_stderr")

IA_2013_stderr <- as.data.frame(read_excel("IA 2013-2014 means stderr.xlsx", sheet = 3, col_types = c("text", rep("numeric", 8)))[c(-1,-272),c(1:5)]) # 188, 252 were missing, added them manually in excel
IA_2013_stderr$PPN <- paste0("SAM", str_pad(as.character(IA_2013_stderr$PPN),3,pad = "0")) # add "SAM" in front of each individual's number
rownames(IA_2013_stderr) <- IA_2013_stderr$PPN
colnames(IA_2013_stderr) <- c("Line", "Palm_stderr", "Stea_stderr", "Oleic_stderr", "Lino_stderr")
IA_2013_stderr[c(4,69,116),2:5] <- NA # remove questionable genotypes

# hist(IA_2013_stderr$Palm_stderr, breaks = 50)
# shapiro.test(IA_2013_stderr$Palm_stderr)

IA_2014_stderr <- as.data.frame(read_excel("IA 2013-2014 means stderr.xlsx", sheet = 4, col_types = c("text", rep("numeric", 8)))[c(-1,-272),c(1:5)]) # 188, 252 were missing, added them manually in excel
IA_2014_stderr$PPN <- paste0("SAM", str_pad(as.character(IA_2014_stderr$PPN),3,pad = "0")) # add "SAM" in front of each individual's number
rownames(IA_2014_stderr) <- IA_2014_stderr$PPN
colnames(IA_2014_stderr) <- c("Line", "Palm_stderr", "Stea_stderr", "Oleic_stderr", "Lino_stderr")
IA_2014_stderr[c(4,69,116),2:5] <- NA # remove questionable genotypes

plot(IA_2014_stderr[-1])
all_stderrs <- rbind(MN_2016_day_1_stderr, MN_2016_day_2_stderr, BC_stderr, GA_stderr, IA_2013_stderr, IA_2014_stderr)

png("fatty_acid_stderrs_relationships.png", width = 1000, height = 1000, res = 100)
#plot(all_stderrs[-1])
chart.Correlation(all_stderrs[-1])
dev.off()

cor.test(all_stderrs$Oleic_stderr, all_stderrs$Lino_stderr)
cor.test(all_stderrs$Oleic_stderr, all_stderrs$Stea_stderr)
cor.test(all_stderrs$Oleic_stderr, all_stderrs$Palm_stderr)


# try GWAS with power transformed stderrs
# IA_2014_stderr_normalized <- transfo(IA_2014_stderr[,2:5])[["Xt"]]
# hist(IA_2014_stderr_normalized[,1], breaks = 50)
# shapiro.test(IA_2014_stderr_normalized[,1])

write.table(MN_2016_day_1_stderr, file = 'MN_2016_day_1_oil_contents_stderr.csv', sep = ",", row.names = F, col.names = T)
write.table(MN_2016_day_2_stderr, file = 'MN_2016_day_2_oil_contents_stderr.csv', sep = ",", row.names = F, col.names = T)
write.table(BC_stderr, file = 'BC_oil_contents_stderr.csv', sep = ",", row.names = F, col.names = T)
write.table(GA_stderr, file = 'GA_oil_contents_stderr.csv', sep = ",", row.names = F, col.names = T)
write.table(IA_2013_stderr, file = 'IA_2013_oil_contents_stderr.csv', sep = ",", row.names = F, col.names = T)
write.table(IA_2014_stderr, file = 'IA_2014_oil_contents_stderr.csv', sep = ",", row.names = F, col.names = T)

#write.table(IA_2014_stderr_normalized, file = 'IA_2014_oil_contents_stderr_normalized.csv', sep = ",", row.names = T, col.names = NA)

IA_2014_stderr_no_FAD <- IA_2014_stderr
IA_2014_stderr_no_FAD[FAD_status$FAD2.1_flag == 1,2:5] <- NA

#write.table(MN_2016_day_1, file = "MN_2016_day_1_oil_contents_no_FAD.csv", sep = ",", row.names = F, col.names = T)
write.table(IA_2014_stderr_no_FAD, file = "IA_2014_oil_contents_stderr_no_FAD.csv", sep = ",", row.names = F, col.names = T)


# all oleic stderrs in one dataframe
oleic_stderrs <- data.frame(BC_stderr = rep(NA, 287))
rownames(oleic_stderrs) <- BC_stderr$Line
oleic_stderrs$BC_stderr <- BC_stderr$Oleic_stderr
oleic_stderrs$MN_2016_day_1_stderr <- MN_2016_day_1_stderr$Oleic_stderr
oleic_stderrs$MN_2016_day_2_stderr <- MN_2016_day_2_stderr$Oleic_stderr
oleic_stderrs$IA_2013_stderr <- IA_2013_stderr$Oleic_stderr
oleic_stderrs$IA_2014_stderr <- IA_2014_stderr$Oleic_stderr
oleic_stderrs$GA_stderr <- GA_stderr$Oleic_stderr

# plot stderrs over different environment
oleic_stderrs$sample <- row.names(oleic_stderrs)



# create parallel coordinate plot of stderrs of oleic
png(filename = "parallel_plot_oleic_stderr.png", width = 2000, height = 1500, res = 120)
ggparcoord(oleic_stderrs[-c(29,186),],
           columns = 1:6, 
           groupColumn = 7, 
           scale = "globalminmax",
           title = "Standard errors of oleic acid",
           showPoints = TRUE, 
           alphaLines = 0.5
) +
theme_bw() +
theme(
    legend.position="none",
    axis.text.x = element_text(angle = 45, vjust = 0.5, size = 12)
)
dev.off()

# create boxplots of stderr of oleic
#boxplot(oleic_stderrs[,1:6])


oleic_stderrs_melt <- melt(oleic_stderrs, na.rm = FALSE)
oleic_stderrs_melt$variable <- as.factor(oleic_stderrs_melt$variable)
oleic_stderrs_melt$sample <- as.factor(oleic_stderrs_melt$sample)


png(filename = "boxplot_oleic_stderr.png", width = 800, height = 800, res = 100)
ggplot(oleic_stderrs_melt, aes(x=variable, y=value)) + 
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 12))+
  xlab("")
dev.off()

# try GxE
library(statgenGxE)
gxe <- statgenSTA::createTD(data = na.omit(oleic_stderrs_melt), genotype = "sample", trial = "variable")

plot(gxe, plotType = "box", traits = "value", orderBy = "ascending")

plot(gxe, plotType = "scatter", traits = "value", orderBy = "ascending")

VarComp <- gxeVarComp(TD = gxe, trait = "value") # with trial as nesting factor: Error: number of levels of each grouping factor must be < number of observations (problems: genotype:trial)
summary(VarComp)

# perform ANOVA on stderr of oleic vs environment. to-do: add GxE

anova <- aov(value~variable + sample, data = oleic_stderrs_melt)
summary(anova) # p-val<2e-16

# biplot of stderrs
BC_stderr$Region <- "BC" # add sample region as a column to easier differentiate them after merging
GA_stderr$Region <- "GA"
MN_2016_day_1_stderr$Region <- "MN_2016_day_1"
MN_2016_day_2_stderr$Region <- "MN_2016_day_2"
IA_2013_stderr$Region <- "IA_2013"
IA_2014_stderr$Region <- "IA_2014"

all_stderrs <- bind_rows(BC_stderr,GA_stderr,MN_2016_day_1_stderr,MN_2016_day_2_stderr, IA_2013_stderr, IA_2014_stderr)
#all_phenotypes <- bind_rows(BC,GA,MN_2016_day_1,MN_2016_day_2, IA_2013, IA_2014)

all_stderrs_melt <- melt(all_stderrs)

all_stderrs_melt$ID <- with(all_stderrs_melt, paste(Region, variable, sep = "_"))

gge <- gge(na.omit(all_stderrs_melt), value ~ Line*ID, scale=T)

png("GGE_stderrs_biplot.png", width = 1000, height = 1000, res = 120)
biplot(gge, origin = 0, hull=F)
dev.off()
png("GGE_stderrs_biplot_addition.png", width = 1000, height = 1000, res = 120)
plot(gge)
dev.off()

biplot3d(gge, origin = 0, hull=F) # pretty fancy
png("GGE_stderrs_biplot_no_genotypes.png", width = 1000, height = 1000, res = 120)
biplot(gge, origin = 0, hull=F, cex.gen = 0)
dev.off()

# create lists with all envs per FA from scratch

linoleic <- data.frame(BC = rep(NA, 287))
rownames(linoleic) <- BC$Line
linoleic$BC <- BC$Lino
linoleic$IA <- IA$Lino
linoleic$GA <- GA$Lino
linoleic$MN_2015 <- MN_2015$Lino
linoleic$MN_2016_day_1 <- MN_2016_day_1$Lino
linoleic$MN_2016_day_2 <- MN_2016_day_2$Lino
linoleic$IA_2013 <- IA_2013$Lino
linoleic$IA_2014 <- IA_2014$Lino

oleic <- data.frame(BC = rep(NA, 287))
rownames(oleic) <- BC$Line
oleic$BC <- BC$Oleic
oleic$IA <- IA$Oleic
oleic$GA <- GA$Oleic
oleic$MN_2015 <- MN_2015$Oleic
oleic$MN_2016_day_1 <- MN_2016_day_1$Oleic
oleic$MN_2016_day_2 <- MN_2016_day_2$Oleic
oleic$IA_2013 <- IA_2013$Oleic
oleic$IA_2014 <- IA_2014$Oleic

stearic <- data.frame(BC = rep(NA, 287))
rownames(stearic) <- BC$Line
stearic$BC <- BC$Stea
stearic$IA <- IA$Stea
stearic$GA <- GA$Stea
stearic$MN_2015 <- MN_2015$Stea
stearic$MN_2016_day_1 <- MN_2016_day_1$Stea
stearic$MN_2016_day_2 <- MN_2016_day_2$Stea
stearic$IA_2013 <- IA_2013$Stea
stearic$IA_2014 <- IA_2014$Stea

palmitic <- data.frame(BC = rep(NA, 287))
rownames(palmitic) <- BC$Line
palmitic$BC <- BC$Palm
palmitic$IA <- IA$Palm
palmitic$GA <- GA$Palm
palmitic$MN_2015 <- MN_2015$Palm
palmitic$MN_2016_day_1 <- MN_2016_day_1$Palm
palmitic$MN_2016_day_2 <- MN_2016_day_2$Palm
palmitic$IA_2013 <- IA_2013$Palm
palmitic$IA_2014 <- IA_2014$Palm

write.table(linoleic, file = "linoleic_all_locations_updated.csv", sep = ",", row.names = T, col.names=NA)
write.table(oleic, file = "oleic_all_locations_updated.csv", sep = ",", row.names = T, col.names = NA)
write.table(stearic, file = "stearic_all_locations_updated.csv", sep = ",", row.names = T, col.names = NA)
write.table(palmitic, file = "palmitic_all_locations_updated.csv", sep = ",", row.names = T, col.names = NA)


# create lists with all envs per FA stderrs from scratch

oleic_stderrs <- data.frame(BC_stderr = rep(NA, 287))
rownames(oleic_stderrs) <- BC_stderr$Line
oleic_stderrs$BC_stderr <- BC_stderr$Oleic_stderr
oleic_stderrs$MN_2016_day_1_stderr <- MN_2016_day_1_stderr$Oleic_stderr
oleic_stderrs$MN_2016_day_2_stderr <- MN_2016_day_2_stderr$Oleic_stderr
oleic_stderrs$IA_2013_stderr <- IA_2013_stderr$Oleic_stderr
oleic_stderrs$IA_2014_stderr <- IA_2014_stderr$Oleic_stderr
oleic_stderrs$GA_stderr <- GA_stderr$Oleic_stderr


linoleic_stderrs <- data.frame(BC_stderr = rep(NA, 287))
rownames(linoleic_stderrs) <- BC_stderr$Line
linoleic_stderrs$BC_stderr <- BC_stderr$Lino_stderr
linoleic_stderrs$MN_2016_day_1_stderr <- MN_2016_day_1_stderr$Lino_stderr
linoleic_stderrs$MN_2016_day_2_stderr <- MN_2016_day_2_stderr$Lino_stderr
linoleic_stderrs$IA_2013_stderr <- IA_2013_stderr$Lino_stderr
linoleic_stderrs$IA_2014_stderr <- IA_2014_stderr$Lino_stderr
linoleic_stderrs$GA_stderr <- GA_stderr$Lino_stderr

stea_stderrs <- data.frame(BC_stderr = rep(NA, 287))
rownames(stea_stderrs) <- BC_stderr$Line
stea_stderrs$BC_stderr <- BC_stderr$Stea_stderr
stea_stderrs$MN_2016_day_1_stderr <- MN_2016_day_1_stderr$Stea_stderr
stea_stderrs$MN_2016_day_2_stderr <- MN_2016_day_2_stderr$Stea_stderr
stea_stderrs$IA_2013_stderr <- IA_2013_stderr$Stea_stderr
stea_stderrs$IA_2014_stderr <- IA_2014_stderr$Stea_stderr
stea_stderrs$GA_stderr <- GA_stderr$Stea_stderr

palmitic_stderrs <- data.frame(BC_stderr = rep(NA, 287))
rownames(palmitic_stderrs) <- BC_stderr$Line
palmitic_stderrs$BC_stderr <- BC_stderr$Palm_stderr
palmitic_stderrs$MN_2016_day_1_stderr <- MN_2016_day_1_stderr$Palm_stderr
palmitic_stderrs$MN_2016_day_2_stderr <- MN_2016_day_2_stderr$Palm_stderr
palmitic_stderrs$IA_2013_stderr <- IA_2013_stderr$Palm_stderr
palmitic_stderrs$IA_2014_stderr <- IA_2014_stderr$Palm_stderr
palmitic_stderrs$GA_stderr <- GA_stderr$Palm_stderr

write.table(linoleic_stderrs, file = "linoleic_stderrs_all_locations.csv", sep = ",", row.names = T, col.names=NA)
write.table(oleic_stderrs, file = "oleic_stderrs_all_locations.csv", sep = ",", row.names = T, col.names = NA)
write.table(stea_stderrs, file = "stearic_stderrs_all_locations.csv", sep = ",", row.names = T, col.names = NA)
write.table(palmitic_stderrs, file = "palmitic_stderrs_all_locations.csv", sep = ",", row.names = T, col.names = NA)

