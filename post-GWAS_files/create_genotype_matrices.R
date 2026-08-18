# create genotype matrices of multivariate significant SNPs

library(vcfR)

# load GEMMA results from post_GWAS_analysis.R

# read VCF to later subset by significant SNPs

vcf <- read.vcfR("../GWAS_inputs/IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf")

# beta
# match df length to the vcf
GWAS_stability_multivariate_temp <- merge(min_p_val_multivariate,GWAS_stability_multivariate, by = c("chr", "ps"), all.x = T, sort = F)


vcf_beta <- vcf[GWAS_stability_multivariate_temp$p_wald.y < (0.05/10937),]

beta_genotype_matrix <- extract.gt(vcf_beta, return.alleles = T)

write.csv(beta_genotype_matrix, "beta_stability_genotype_matrix.csv")

# multivariate alpha, using the minimum p-values across environments with the less strict threshold
vcf_alpha <- vcf[min_p_val_multivariate_stderr$p_wald < (0.05/10937),]

alpha_genotype_matrix <- extract.gt(vcf_alpha, return.alleles = T)

write.csv(alpha_genotype_matrix, "alpha_stability_genotype_matrix.csv")

# multivariate mean FA contents, no FAD omission, using the minimum p-values across environments with the less strict threshold
vcf_FA_mean <- vcf[min_p_val_multivariate$p_wald < (0.05/10937),]

FA_mean_genotype_matrix <- extract.gt(vcf_FA_mean, return.alleles = T)

write.csv(FA_mean_genotype_matrix, "FA_mean_stability_genotype_matrix.csv")

# multivariate mean FA contents, FAD mutant omission, using the minimum p-values across environments with the less strict threshold
# match df length to the vcf
min_p_val_multivariate_no_FAD_temp <- merge(min_p_val_multivariate,min_p_val_multivariate_no_FAD, by = c("chr", "ps"), all.x = T, sort = F)

vcf_FA_mean_no_FAD <- vcf[min_p_val_multivariate_no_FAD_temp$p_wald.y < (0.05/10937),]

FA_mean_no_FAD_genotype_matrix <- extract.gt(vcf_FA_mean_no_FAD, return.alleles = T)

write.csv(FA_mean_no_FAD_genotype_matrix, "FA_mean_no_FAD_stability_genotype_matrix.csv")
