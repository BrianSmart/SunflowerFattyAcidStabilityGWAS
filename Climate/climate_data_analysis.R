# RDA plot of fatty acids and climate variables
# create 4 different RDA plots for each fatty acid

library(tidyverse)
library(ggplot2)
library(corrplot)
library(vegan)
library(ggfortify)
library(ggpubr)
#library(ggvegan)
library(tidyr)
library(dplyr)
library(purrr)
library(psych)

setwd("~/projects/GWAS_Markus/Climate")

climates <- read.csv("climate_data.csv", row.names = 1)
climates <- climates[order(rownames(climates)), ] # make sure both inputs match
climates <- scale(climates) # scale climate input
res <- cor(climates) # correlations between climate variables
corrplot(res, type = "upper", order = "alphabet", 
         tl.col = "black", tl.srt = 45)

test <- prcomp(climates)
autoplot(test, label=T, label.repel=T)
test$rotation
biplot(test)

# load phenotype data and perform RDA
# turns out there are too many dimensions for this, common garden stuff only works because every genotype corresponds to an environment
# 
# oleic <- t(na.omit(read.csv("../Phenotype_Analysis/oleic_all_locations_updated.csv", row.names = 1)[-c(2,4,7)])) # omit IA_2013 until I get the climate data and IA_2010 for missing data
# 
# rownames(oleic) <- c("BC", "GA", "MN_2016_early", "MN_2016_early", "IA_2014")
# oleic <- oleic[order(rownames(oleic)), ] # make sure both inputs match
# 
# rda <- rda(oleic, climates[-c(3,5),], scale=F)
# plot(rda,type="text", scaling=0)
# 
# 
# linoleic <- t(na.omit(read.csv("../Phenotype_Analysis/linoleic_all_locations_updated.csv", row.names = 1)[-c(2,4,7)])) # omit IA_2013 until I get the climate data and IA_2010 for missing data
# 
# rownames(linoleic) <- c("BC", "GA", "MN_2016_early", "MN_2016_early", "IA_2014")
# linoleic <- linoleic[order(rownames(linoleic)), ] # make sure both inputs match
# 
# rda <- rda(linoleic, climates[-c(3,5),], scale=F)
# plot(rda,type="text", scaling=0)
# autoplot(rda)
# 
# palm <- t(na.omit(read.csv("../Phenotype_Analysis/palmitic_all_locations_updated.csv", row.names = 1)[-c(2,4,7)])) # omit IA_2013 until I get the climate data and IA_2010 for missing data
# 
# rownames(palm) <- c("BC", "GA", "MN_2016_early", "MN_2016_early", "IA_2014")
# palm <- palm[order(rownames(palm)), ] # make sure both inputs match
# 
# rda <- rda(palm, climates[-c(3,5),], scale=F)
# plot(rda,type="text", scaling=0)
# 
# stearic <- t(na.omit(read.csv("../Phenotype_Analysis/stearic_all_locations_updated.csv", row.names = 1)[-c(2,4,7)])) # omit IA_2013 until I get the climate data and IA_2010 for missing data
# 
# rownames(stearic) <- c("BC", "GA", "MN_2016_early", "MN_2016_early", "IA_2014")
# stearic <- stearic[order(rownames(stearic)), ] # make sure both inputs match
# 
# rda <- rda(stearic, climates[-c(3,5),], scale=F)
# plot(rda,type="text", scaling=0)
#
# try all four phenotypes at once => doesn't work this way
# test <- tibble(oleic = oleic, linoleic = linoleic, stearic = stearic, palm = palm)
# rda <- rda(test, climates[-c(3,5),], scale=F)
# plot(rda,type="text", scaling=0)

# look for correlations between mean fatty acids and stderr and climate variables
# only use samples present in all environments to prevent bias
oleic <- read.csv("../Phenotype_Analysis/oleic_all_locations_updated.csv", row.names = 1)[-c(2)] # omit IA_2010 for missing data
linoleic <- read.csv("../Phenotype_Analysis/linoleic_all_locations_updated.csv", row.names = 1)[-c(2)] # omit IA_2010 for missing data
palm <- read.csv("../Phenotype_Analysis/palmitic_all_locations_updated.csv", row.names = 1)[-c(2)] # omit IA_2010 for missing data
stearic <- read.csv("../Phenotype_Analysis/stearic_all_locations_updated.csv", row.names = 1)[-c(2)] # omit IA_2010 for missing data
FA_means <- data.frame(oleic = apply(t(na.omit(oleic)), 1, mean), linoleic = apply(t(na.omit(linoleic)), 1, mean), palm = apply(t(na.omit(palm)), 1, mean), stearic = apply(t(na.omit(stearic)), 1, mean))
FA_means <- FA_means[order(rownames(FA_means)), ] # make sure both inputs match

# split by FAD mutation status
FAD_mutation <- read.csv("../Phenotype_Analysis/FAD2-1_mutant_status_renamed.csv")

oleic_no_FAD <- oleic[FAD_mutation$FAD2.1_flag == F,]
linoleic_no_FAD <- linoleic[FAD_mutation$FAD2.1_flag == F,]
palm_no_FAD <- palm[FAD_mutation$FAD2.1_flag == F,]
stearic_no_FAD <- stearic[FAD_mutation$FAD2.1_flag == F,]
FA_means_no_FAD <- data.frame(oleic = apply(t(na.omit(oleic_no_FAD)), 1, mean), linoleic = apply(t(na.omit(linoleic_no_FAD)), 1, mean), palm = apply(t(na.omit(palm_no_FAD)), 1, mean), stearic = apply(t(na.omit(stearic_no_FAD)), 1, mean))
FA_means_no_FAD <- FA_means_no_FAD[order(rownames(FA_means_no_FAD)), ] # make sure both inputs match

# only one FAD mutant is present in all environments 
oleic_FAD <- oleic[FAD_mutation$FAD2.1_flag == T,]
linoleic_FAD <- linoleic[FAD_mutation$FAD2.1_flag == T,]
palm_FAD <- palm[FAD_mutation$FAD2.1_flag == T,]
stearic_FAD <- stearic[FAD_mutation$FAD2.1_flag == T,]
FA_means_FAD <- data.frame(oleic = apply(t(na.omit(oleic_FAD)), 1, mean), linoleic = apply(t(na.omit(linoleic_FAD)), 1, mean), palm = apply(t(na.omit(palm_FAD)), 1, mean), stearic = apply(t(na.omit(stearic_FAD)), 1, mean))
FA_means_FAD <- FA_means_FAD[order(rownames(FA_means_FAD)), ] # make sure both inputs match


cor_p_vals <- corr.test(FA_means_no_FAD,climates[,order(colnames(climates))][-3,]) # get p-values for each correlation
colors <- ifelse(c(cor_p_vals$p < 0.05),"red", "black") # mark significant correlations red

#rownames(FA_means) <- c("BC", "IA_2010", "GA", "MN_2015", "MN_2016_day_1", "MN_2016_day_2", "IA_2014")
cor <- cor(FA_means_no_FAD,climates[,order(colnames(climates))][-3,]) # correlations are computed, climates are ordered by alphabet

png("correlation_plot_no_FAD.png", width = 2000, height = 400)
corrplot(cor, insig = "pch", method = "color", addCoef.col=colors, 
         tl.col = "black", tl.srt = 45,tl.cex = 1/par("cex"),
         cl.cex = 1/par("cex"), addCoefasPercent = TRUE)
dev.off()

# to-do: with stderrs
oleic_stderrs <- t(na.omit(read.csv("../Phenotype_Analysis/oleic_stderrs_all_locations.csv", row.names = 1)[-c(4)])) # omit IA_2013 until I get the climate data
linoleic_stderrs <- t(na.omit(read.csv("../Phenotype_Analysis/linoleic_stderrs_all_locations.csv", row.names = 1)[-c(4)])) # omit IA_2013 until I get the climate data
palm_stderrs <- t(na.omit(read.csv("../Phenotype_Analysis/palmitic_stderrs_all_locations.csv", row.names = 1)[-c(4)])) # omit IA_2013 until I get the climate data
stearic_stderrs <- t(na.omit(read.csv("../Phenotype_Analysis/stearic_stderrs_all_locations.csv", row.names = 1)[-c(4)])) # omit IA_2013 until I get the climate data
FA_stderrs_means <- data.frame(oleic_stderrs = apply(oleic_stderrs, 1, mean), linoleic_stderrs = apply(linoleic_stderrs, 1, mean), palm_stderrs = apply(palm_stderrs, 1, mean), stearic_stderrs = apply(stearic_stderrs, 1, mean))
FA_stderrs_means <- FA_stderrs_means[order(rownames(FA_stderrs_means)), ] # make sure both inputs match

cor_p_vals <- corr.test(FA_stderrs_means ,climates[,order(colnames(climates))][-c(3,5),]) # get p-values for each correlation
colors <- ifelse(c(cor_p_vals$p < 0.05),"red", "black") # mark significant correlations red

cor <- cor(FA_stderrs_means,climates[,order(colnames(climates))][-c(3,5),]) # correlations are computed, climates are ordered by alphabet

png("correlation_stderrs_plot.png", width = 2000, height = 400)
corrplot(cor, insig = "pch", method = "color", addCoef.col=colors, 
         tl.col = "black", tl.srt = 45,tl.cex = 1/par("cex"),
         cl.cex = 1/par("cex"), addCoefasPercent = TRUE)
dev.off()

