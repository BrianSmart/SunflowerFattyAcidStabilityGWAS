# plot gene density to maybe see where the centromere lies

setwd("~/Genotype_Analysis")

library(ape)
library(zoo)
library(ggplot2)
library(ggpubr)
library(dplyr)

gff <- read.gff("../Post-GWAS_analysis/HAN412_Eugene_curated_v1_1.gff3")
genes <- subset(gff, type == "gene")

ggplot(subset(genes, seqid == "Ha412HOChr17"), aes(start)) +
  geom_histogram(fill= "blue",alpha=0.5, bins = 80) +
  theme_bw() +
  #ggtitle("Chr05") + 
  theme(text = element_text(size = 20))