# analysis of the Chr 5 introgression identified by Todesco et al.
library(readxl)
#library(GMMAT)
library(scales)
library(vcfR)
library(stringr)
library(snpStats)
setwd("~/projects/GWAS_FattyAcid_Stability/Introgression_Analysis/")


chr5introgressions <- read_xlsx("ann05.01 haplotypes SAM.xlsx", col_names = F)[,c(1,6)]
# GA_stderrs <- read.csv("../Phenotype_Analysis/GA_oil_contents_stderr.csv", header = T)
# MN_2016_early_stderrs <- read.csv("../Phenotype_Analysis/MN_2016_day_1_oil_contents_stderr.csv", header = T)

#GEMMA_kinship <- read.csv("../GWAS_inputs/SAM_Kinship_new.cXX.txt", sep = "\t", header = F)
#introgression_GA_stderrs <- merge(chr5introgressions, GA_stderrs, by.x = "...1", by.y = "Line")
#introgression_GA_stderrs$introgression_scaled <- rescale(introgression_GA_stderrs$...6)

# didn't get how to add the kinship matrix, just using 2 PCs instead
PCs <- read.csv("../PCA.R/PC_covariates_GEMMA.csv", header = F)[-1]
introgression_GA_stderrs$PC1 <- PCs$V2
introgression_GA_stderrs$PC2 <- PCs$V3

# introgression presence vs. oleic stderr

lm <- aov(Stea_stderr ~ ...6 + PC1 + PC2, data = na.omit(introgression_GA_stderrs))
summary(lm)
lm$coefficients

plot(introgression_GA_stderrs$...6, introgression_GA_stderrs$Stea_stderr, xlab = "introgression quantity", ylab = "GA stearic stderrs")
abline(lm)

# check if significant SNPs are in LD with the Chr5 introgression
GWAS_GA_multivariate_stderr <- read.delim("../GWAS_outputs/GEMMA_GA_stderr.assoc.txt", header = T)

vcf <- read.vcfR("../GWAS_inputs/IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf")

vcf_GA_stderr <- vcf[GWAS_GA_multivariate_stderr$p_wald < (0.05/2208147),]

GA_stderr_genotype_matrix <- extract.gt(vcf_GA_stderr, return.alleles = F, as.numeric = F) # as numeric = T only returns 0 or 1
dosages <- sapply(GA_stderr_genotype_matrix, function(x) sum(as.numeric(unlist(strsplit(x,split = "|"))[-2])))
dosage_matrix <- matrix(dosages, ncol = 287)
dosage_matrix_introgression <- rbind(dosage_matrix, chr5introgressions$...6)
rownames(dosage_matrix_introgression) <- c(rownames(GA_stderr_genotype_matrix), "chr5_introgression")
colnames(dosage_matrix_introgression) <- colnames(GA_stderr_genotype_matrix)

geno_matrix <-new("SnpMatrix", t(dosage_matrix_introgression))
ld_all <- ld(geno_matrix[,1:328], stats = "R.squared", depth = 329)
ld_snps <- ld(geno_matrix[,11:201], stats = "R.squared", depth = 329)
ld_chr5 <- ld(geno_matrix[,11:201], geno_matrix[,329], stats = "R.squared", depth = 329)

png("rsq_all_significant_SNPs_GA_stderrs.png", res=120, width = 1500, height = 1500)
image(ld_all, lwd=0, cuts=9, colorkey=T, col = rev(heat.colors(20)), main="rsq all SNPs")
dev.off()

png("rsq_chr_05_significant_SNPs_GA_stderrs.png", res=120, width = 1500, height = 1500)
image(ld_snps, lwd=0, cuts=9, colorkey=TRUE, col = rev(heat.colors(20)), main="rsq Chr 5")
dev.off()

png("rsq_chr_05_significant_SNPs_GA_stderrs_vs_introgression.png", res=120, width = 1000, height = 500)
image(ld_chr5, lwd=0, colorkey=TRUE, col = rev(heat.colors(20)), main="rsq Chr 5 SNPs vs introgression")
dev.off()
write.csv(ld_chr5, file = "rsq_introgression.csv")
