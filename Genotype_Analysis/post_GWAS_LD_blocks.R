library(VariantAnnotation)
library(updog)
library(ldsep)


setwd("D:/Sunflower/New_Plots")
# subset for significant regions
GWAS_stability_multivariate <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_stability.assoc.txt", header = T)
GWAS_stability_multivariate <- GWAS_stability_multivariate[GWAS_stability_multivariate$p_wald<0.05/length(GWAS_stability_multivariate$chr),]


setwd("D:/Sunflower_stability_GWAS/Genotype_Analysis/")

# subset as in https://www.biostars.org/p/170965/

#test <- vcfR::read.vcfR("D:/Sunflower_stability_GWAS/Genotype_Analysis/IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf", nrows = 100)
test <- vcfR::read.vcfR("variants_of_interest.vcf", nrows = 100)
# rngs <- GRanges("1", IRanges(c(4426097, 214822553)))
# param <- ScanVcfParam(which=rngs) 

VCF <- readVcf("variants_of_interest.vcf")



sizemat <- geno(VCF)$GS
refmat <- geno(VCF)$RA

ploidy <- 2
mout <- multidog(refmat = refmat, 
                 sizemat = sizemat, 
                 ploidy = ploidy, 
                 model = "norm")