# analysis of Introgression GEMMA output
library(qqman)
library(GWASTools)
library(SeqArray)
library(ape)
library(zoo)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(qqplotr)
library(stringr)
library(naniar)
library(readxl)
library(arsenal)
library(CMplot)
library(tidyverse)

setwd("~/projects/GWAS_FattyAcid_Stability/Introgression_Analysis/")

# look at a representative/interesting subset of results

GWAS_IA_2014_introgression <- read.delim("GWAS_outputs/GEMMA_IA_2014_Introgressions.assoc.txt", header = T)
GWAS_IA_2014_no_FAD_introgression <- read.delim("GWAS_outputs/GEMMA_IA_2014_no_FAD_Introgressions.assoc.txt", header = T)
GWAS_IA_2014_stderr_introgression <- read.delim("GWAS_outputs/GEMMA_IA_2014_stderr_Introgressions.assoc.txt", header = T)
GWAS_MN_2016_early_introgression <- read.delim("GWAS_outputs/GEMMA_MN_2016_early_Introgressions.assoc.txt", header = T)
GWAS_MN_2016_early_no_FAD_introgression <- read.delim("GWAS_outputs/GEMMA_MN_2016_early_no_FAD_Introgressions.assoc.txt", header = T)
GWAS_MN_2016_early_stderr_introgression <- read.delim("GWAS_outputs/GEMMA_MN_2016_early_stderr_Introgressions.assoc.txt", header = T)
GWAS_GA_stderr_introgression <- read.delim("GWAS_outputs/GEMMA_GA_stderr_Introgressions.assoc.txt", header = T)

png("Manhattan_plot_Introgressions.png", res=120, width = 1000, height = 3500)

layout(matrix(c(1,2,3,4,5,6,7), byrow = F),
       heights =c(6,6,6,6,6,6,7))

manhattan(GWAS_IA_2014_introgression, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA 2014"),
          suggestiveline = FALSE, genomewideline = -log10(0.05/266), )
manhattan(GWAS_IA_2014_no_FAD_introgression, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA 2014 no FAD"),
          suggestiveline = FALSE, genomewideline = -log10(0.05/266))
manhattan(GWAS_IA_2014_stderr_introgression, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA 2014 stderr"),
          suggestiveline = FALSE, genomewideline = -log10(0.05/266))
manhattan(GWAS_MN_2016_early_introgression, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_early"),
          suggestiveline = FALSE, genomewideline = -log10(0.05/266), )
manhattan(GWAS_MN_2016_early_no_FAD_introgression, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_early no FAD"),
          suggestiveline = FALSE, genomewideline = -log10(0.05/266))
manhattan(GWAS_MN_2016_early_stderr_introgression, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_early stderr"),
          suggestiveline = FALSE, genomewideline = -log10(0.05/266))
manhattan(GWAS_GA_stderr_introgression, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("GA stderr"),
          suggestiveline = FALSE, genomewideline = -log10(0.05/266))
dev.off()
