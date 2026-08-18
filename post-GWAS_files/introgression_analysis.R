# compare introgressions found by Todesco et al. 2020 with significant SNPs from multivariate GWAS
library(ggplot2)
library(ggpubr)
library(dplyr)
library(qqplotr)
library(stringr)
library(naniar)
library(readxl)
library(tidyverse)
library(xtable)
library(GenomicRanges)
library(Repitools)
library(GenomicDistributions)

setwd("~/projects/GWAS_FattyAcid_Stability/post-GWAS_files")

# load and subset introgression list
introgressions <- read.csv("SAM_introgression.pcadmix.txt", sep = "\t")
introgressions <- introgressions[order(introgressions$CHR),]
introgressions$CHR <- as.integer(gsub("Ha412HOChr", "", introgressions$CHR))

introgressions_freq <- plyr::count(introgressions[,c(2,3,5)], c("start", "end", "CHR")) # count how often each introgression is present
introgressions_freq$width <- introgressions_freq$end - introgressions_freq$start
introgressions_freq_gt_10 <- introgressions_freq[introgressions_freq$freq > 28,] # MAF > 0.05, as A and B allles of the same samples were counted double

write.csv(introgressions_freq_gt_10, file = "SAM_introgression.pcadmix_gt_14.csv")


# how much of the genome is covered by introgressions?
# extract the chromosome length from the fasta
sunflower_chromosomes <- getChromSizesFromFasta("Ha412HOv2.0-20181130.fasta.gz")[1:17]
names(sunflower_chromosomes) <- 1:17

# create genomic ranges object from introgressions
introgressions_GR <- makeGRangesFromDataFrame(introgressions_freq, start.field = "start", end.field = "end")
introgression_coverage <- coverage(introgressions_GR) # extract coverages from GR object, resulting in coverage per base
sum(introgression_coverage>0)/sunflower_chromosomes*100 # percentage of basepairs per chromosome covered with introgressions


# load GWAS results, subset significant SNPs
min_p_val_multivariate <- read.csv("min_p_val_multivariate.csv")
min_p_val_multivariate_subset <- min_p_val_multivariate[min_p_val_multivariate$p_wald < 0.05/length(min_p_val_multivariate$chr),] # subset for significant SNPs. Stricter threshold due to p-val inflation
#min_p_val_multivariate_subset$end <- min_p_val_multivariate_subset$ps

min_p_val_multivariate_no_FAD <- read.csv("min_p_val_multivariate_no_FAD.csv")
min_p_val_multivariate_no_FAD_subset <- min_p_val_multivariate_no_FAD[min_p_val_multivariate_no_FAD$p_wald < 0.05/length(min_p_val_multivariate_no_FAD$chr),] # subset for significant SNPs. Stricter threshold for consistency

min_p_val_multivariate_stderr <- read.csv("min_p_val_multivariate_stderr.csv")
min_p_val_multivariate_stderr_subset <- min_p_val_multivariate_stderr[min_p_val_multivariate_stderr$p_wald < 0.05/length(min_p_val_multivariate_stderr$chr),] # subset for significant SNPs. Stricter threshold for consistency

GWAS_stability_multivariate <- read.delim("../GWAS_outputs/GEMMA_stability.assoc.txt", header = T)
stability_p_val_subset <- GWAS_stability_multivariate[GWAS_stability_multivariate$p_wald < 0.05/length(GWAS_stability_multivariate$chr),]

# # use the GenomicRanges package to compare both datasets
# min_p_val_multivariate_GR <- makeGRangesFromDataFrame(na.omit(min_p_val_multivariate_subset), start.field = "ps", end.field = "ps", keep.extra.columns = T)
# min_p_val_multivariate_no_FAD_GR <- makeGRangesFromDataFrame(na.omit(min_p_val_multivariate_no_FAD_subset), start.field = "ps", end.field = "ps", keep.extra.columns = T)
# min_p_val_multivariate_stderr_GR <- makeGRangesFromDataFrame(na.omit(min_p_val_multivariate_stderr_subset), start.field = "ps", end.field = "ps", keep.extra.columns = T)
# stability_p_val_GR <- makeGRangesFromDataFrame(na.omit(stability_p_val_subset), start.field = "ps", end.field = "ps", keep.extra.columns = T)
# 
#
# # SNPs in introgressions
# overlap_all <- subsetByOverlaps(introgressions_GR,min_p_val_multivariate_GR, type = "any")
# test <- pintersect(overlap_all)
# overlap_introgressions_all <- annoGR2DF(overlap_all)
# 
# overlap_no_FAD <- subsetByOverlaps(introgressions_GR,min_p_val_multivariate_no_FAD_GR)
# overlap_no_FAD_introgressions_all <- annoGR2DF(overlap_no_FAD)
# 
# overlap_stderr <- subsetByOverlaps(introgressions_GR,min_p_val_multivariate_stderr_GR)
# overlap_stderr_introgressions_all <- annoGR2DF(overlap_stderr)


# try another way of getting overlaps without losing any columns
merged_df <- merge(introgressions_freq_gt_10, min_p_val_multivariate_subset, by=NULL) # cross-join, results in each row configuration being output
overlap_all <- subset(merged_df, chr == CHR & ps > start & ps < end)[-6] # subset for significant SNPs within introgressions
overlap_all_count <- plyr::count(overlap_all[-4], vars=c("start", "end", "CHR"))[,c(3,1,2,4)] # count how many SNPs are in each window, reorder rows
overlap_all_count <- overlap_all_count[order(overlap_all_count$CHR, overlap_all_count$start),] # sort by position
overlap_all_count$length <- overlap_all_count$end - overlap_all_count$start # add window length

merged_df <- merge(introgressions_freq_gt_10, min_p_val_multivariate_no_FAD_subset, by=NULL) # cross-join, results in each row configuration being output
overlap_no_FAD <- subset(merged_df, chr == CHR & ps > start & ps < end)[-6] # subset for significant SNPs within introgressions
overlap_no_FAD_count <- plyr::count(overlap_no_FAD[-4], vars=c("start", "end", "CHR"))[,c(3,1,2,4)] # count how many SNPs are in each window, reorder rows
overlap_no_FAD_count <- overlap_no_FAD_count[order(overlap_no_FAD_count$CHR, overlap_no_FAD_count$start),] # sort by position
overlap_no_FAD_count$length <- overlap_no_FAD_count$end - overlap_no_FAD_count$start # add window length

merged_df <- merge(introgressions_freq_gt_10, min_p_val_multivariate_stderr_subset, by=NULL) # cross-join, results in each row configuration being output
overlap_stderr <- subset(merged_df, chr == CHR & ps > start & ps < end)[-6] # subset for significant SNPs within introgressions
overlap_stderr_count <- plyr::count(overlap_stderr[-4], vars=c("start", "end", "CHR"))[,c(3,1,2,4)] # count how many SNPs are in each window, reorder rows
overlap_stderr_count <- overlap_stderr_count[order(overlap_stderr_count$CHR, overlap_stderr_count$start),] # sort by position
overlap_stderr_count$length <- overlap_stderr_count$end - overlap_stderr_count$start # add window length

merged_df <- merge(introgressions_freq_gt_10, stability_p_val_subset, by=NULL) # cross-join, results in each row configuration being output
overlap_beta <- subset(merged_df, chr == CHR & ps > start & ps < end)[-6] # subset for significant SNPs within introgressions
overlap_beta_count <- plyr::count(overlap_beta[-4], vars=c("start", "end", "CHR"))[,c(3,1,2,4)] # count how many SNPs are in each window, reorder rows
overlap_beta_count <- overlap_beta_count[order(overlap_beta_count$CHR, overlap_beta_count$start),] # sort by position
overlap_beta_count$length <- overlap_beta_count$end - overlap_beta_count$start # add window length

# write results to csv
write.csv(overlap_all_count, file = "introgression_SNP_overlap_all_samples.csv")
write.csv(overlap_no_FAD_count, file = "introgression_SNP_overlap_no_FAD.csv")
write.csv(overlap_stderr_count, file = "introgression_SNP_overlap_stderr.csv")
write.csv(overlap_beta_count, file = "introgression_SNP_overlap_beta.csv")
write.csv(introgressions_freq[,c(3,1,2,5,4)], file = "SAM_introgression.pcadmix_simplified.csv")

# try plotting introgressions and SNPs on the genome using GenomicDistributions package
sunflower_chromosomes <- getChromSizesFromFasta("Ha412HOv2.0-20181130.fasta.gz")[1:17]
names(sunflower_chromosomes) <- 1:17
#sunflower_chr_df <- data.frame(chr=1:17, start = 1, end = sunflower_chromosomes)
bins <- getGenomeBins(sunflower_chromosomes,binCount = 100000)

queryList = GRangesList(introgressions = introgressions_GR,SNPs = min_p_val_multivariate_GR, compress = F)

x=calcChromBins(queryList, bins)
plotChromBins(x)
