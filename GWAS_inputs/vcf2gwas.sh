#!/bin/bash

# !!! At this point, vcf2gwas cannot perform multivariate analyses. Instead, it will only use the first phenotype for the analysis !!!
# vcf2gwas used on all eight environment phenotypes. The results are always stored in a folder called output. 
# To keep results in unique folders, the output folder is renamed after every run, overwriting results folders of the same name!

## -lmm and -gk isn't allowed together. Kinship computed with GEMMA beforehand
##vcf2gwas -v ../GATK_SAM_population/IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf -pf BC_oil_contents.csv -k SAM_Kinship_new.cXX.txt -ap -lmm -m -T 36

# BC phenotype
#vcf2gwas -v ../GATK_SAM_population/IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf -pf BC_oil_contents.csv -k SAM_Kinship_new.cXX.txt -ap -lmm -T 36
#mv -f Output/ vcf2gwas_BC

# GA phenotype
#vcf2gwas -v ../GATK_SAM_population/IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf -pf GA_oil_contents.csv -k SAM_Kinship_new.cXX.txt -ap -lmm -T 36
#mv -f Output/ vcf2gwas_GA

# IA phenotype. No significant SNPs
#vcf2gwas -v ../GATK_SAM_population/IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf -pf IA_oil_contents.csv -k SAM_Kinship_new.cXX.txt -ap -lmm -T 36
#mv -f Output/ vcf2gwas_IA

# MN_2015 phenotype
#vcf2gwas -v ../GATK_SAM_population/IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf -pf MN_2015_oil_contents.csv -k SAM_Kinship_new.cXX.txt -ap -lmm -T 36
#mv -f Output/ vcf2gwas_MN_2015

# MN_2016_day_1 phenotype
#vcf2gwas -v ../GATK_SAM_population/IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf -pf MN_2016_day_1_oil_contents.csv -k SAM_Kinship_new.cXX.txt -ap -lmm -T 36
#mv -f Output/ vcf2gwas_MN_2016_day_1

# MN_2016_day_2 phenotype
#vcf2gwas -v ../GATK_SAM_population/IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf -pf MN_2016_day_2_oil_contents.csv -k SAM_Kinship_new.cXX.txt -ap -lmm -T 36
#mv -f Output/ vcf2gwas_MN_2016_day_2

# IA_2013 phenotype
#vcf2gwas -v ../GATK_SAM_population/IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf -pf IA_2013_oil_contents.csv -k SAM_Kinship_new.cXX.txt -ap -lmm -T 36
#mv -f Output/ vcf2gwas_IA_2013

# IA_2014 phenotype
#vcf2gwas -v ../GATK_SAM_population/IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf -pf IA_2014_oil_contents.csv -k SAM_Kinship_new.cXX.txt -ap -lmm -T 36
#mv -f Output/ vcf2gwas_IA_2014


# results from omitting mutants and using a covariate are nearly identical, with omission p-vals tend to be slightly higher
# Add presence of FAD2-mutation as a covariate. On IA_2014
#vcf2gwas -v ../GATK_SAM_population/IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf -pf IA_2014_oil_contents.csv -k SAM_Kinship_new.cXX.txt -ap -cf FAD2-1_mutant_status_renamed.csv -c 1 -lmm -T 36
#mv -f Output/ vcf2gwas_IA_2014_FAD_covariate

# Remove FAD2-mutants from analysis entirely. On IA_2014
#vcf2gwas -v ../GATK_SAM_population/IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf -pf IA_2014_oil_contents_no_FAD.csv -k SAM_Kinship_new.cXX.txt -ap -lmm -T 36
#mv -f Output/ vcf2gwas_IA_2014_no_FAD_mutants

# Add presence of FAD2-mutation as a covariate. On GA
#vcf2gwas -v ../GATK_SAM_population/IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf -pf GA_oil_contents.csv -k SAM_Kinship_new.cXX.txt -ap -cf FAD2-1_mutant_status_renamed.csv -c 1 -lmm -T 36
#mv -f Output/ vcf2gwas_GA_FAD_covariate

# Add presence of FAD2-mutation as a covariate. On BC
#vcf2gwas -v ../GATK_SAM_population/IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf -pf BC_oil_contents.csv -k SAM_Kinship_new.cXX.txt -ap -cf FAD2-1_mutant_status_renamed.csv -c 1 -lmm -T 36
#mv -f Output/ vcf2gwas_BC_FAD_covariate

# Add presence of FAD2-mutation as a covariate. On IA
#vcf2gwas -v ../GATK_SAM_population/IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf -pf IA_oil_contents.csv -k SAM_Kinship_new.cXX.txt -ap -cf FAD2-1_mutant_status_renamed.csv -c 1 -lmm -T 36
#mv -f Output/ vcf2gwas_IA_FAD_covariate

# Add presence of FAD2-mutation as a covariate. On MN_2015
#vcf2gwas -v ../GATK_SAM_population/IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf -pf MN_2015_oil_contents.csv -k SAM_Kinship_new.cXX.txt -ap -cf FAD2-1_mutant_status_renamed.csv -c 1 -lmm -T 36
#mv -f Output/ vcf2gwas_MN_2015_FAD_covariate

# Add presence of FAD2-mutation as a covariate. On MN_2016_day_1
#vcf2gwas -v ../GATK_SAM_population/IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf -pf MN_2016_day_1_oil_contents.csv -k SAM_Kinship_new.cXX.txt -ap -cf FAD2-1_mutant_status_renamed.csv -c 1 -lmm -T 36
#mv -f Output/ vcf2gwas_MN_2016_day_1_FAD_covariate

# Add presence of FAD2-mutation as a covariate. On MN_2016_day_2
#vcf2gwas -v ../GATK_SAM_population/IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf -pf MN_2016_day_2_oil_contents.csv -k SAM_Kinship_new.cXX.txt -ap -cf FAD2-1_mutant_status_renamed.csv -c 1 -lmm -T 36
#mv -f Output/ vcf2gwas_MN_2016_day_2_FAD_covariate

# Add presence of FAD2-mutation as a covariate. On IA_2013
#vcf2gwas -v ../GATK_SAM_population/IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf -pf IA_2013_oil_contents.csv -k SAM_Kinship_new.cXX.txt -ap -cf FAD2-1_mutant_status_renamed.csv -c 1 -lmm -T 36
#mv -f Output/ vcf2gwas_IA_2013_FAD_covariate



# BC standard error/stability phenotype
#vcf2gwas -v ../GATK_SAM_population/IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf -pf BC_oil_contents_stderr.csv -k SAM_Kinship_new.cXX.txt -ap -lmm -T 36
#mv -f Output/ vcf2gwas_BC_stderr

# GA standard error/stability phenotype
#vcf2gwas -v ../GATK_SAM_population/IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf -pf GA_oil_contents_stderr.csv -k SAM_Kinship_new.cXX.txt -ap -lmm -T 36
#mv -f Output/ vcf2gwas_GA_stderr

# MN_2016_day_1 standard error/stability phenotype
#vcf2gwas -v ../GATK_SAM_population/IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf -pf MN_2016_day_1_oil_contents_stderr.csv -k SAM_Kinship_new.cXX.txt -ap -lmm -T 36
#mv -f Output/ vcf2gwas_MN_2016_day_1_stderr

# MN_2016_day_2 standard error/stability phenotype
#vcf2gwas -v ../GATK_SAM_population/IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf -pf MN_2016_day_2_oil_contents_stderr.csv -k SAM_Kinship_new.cXX.txt -ap -lmm -T 36
#mv -f Output/ vcf2gwas_MN_2016_day_2_stderr

# IA_2013 standard error/stability phenotype
#vcf2gwas -v ../GATK_SAM_population/IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf -pf IA_2013_oil_contents_stderr.csv -k SAM_Kinship_new.cXX.txt -ap -lmm -T 36
#mv -f Output/ vcf2gwas_IA_2013_stderr

# IA_2014 standard error/stability phenotype
#vcf2gwas -v ../GATK_SAM_population/IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf -pf IA_2014_oil_contents_stderr.csv -k SAM_Kinship_new.cXX.txt -ap -lmm -T 36
#mv -f Output/ vcf2gwas_IA_2014_stderr

# IA_2014 standard error/stability phenotype without FAD-mutants
#vcf2gwas -v ../GATK_SAM_population/IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf -pf IA_2014_oil_contents_stderr_no_FAD.csv -k SAM_Kinship_new.cXX.txt -ap -lmm -T 36
#mv -f Output/ vcf2gwas_IA_2014_stderr_no_FAD

## IA_2014 standard error/stability phenotype normalized
##vcf2gwas -v ../GATK_SAM_population/IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf -pf IA_2014_oil_contents_stderr_normalized.csv -k SAM_Kinship_new.cXX.txt -ap -lmm -T 36
##mv -f Output/ vcf2gwas_IA_2014_stderr_norm

# phenotype stability over all environments as input
#vcf2gwas -v ../GATK_SAM_population/IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf -pf oil_stabilities_normalized.csv -k SAM_Kinship_new.cXX.txt -ap -lmm -T 36
#mv -f Output/ vcf2gwas_oil_stability

# phenotype stderr stability over all environments as input
#vcf2gwas -v ../GATK_SAM_population/IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf -pf oil_stderrs_stabilities_normalized.csv -k SAM_Kinship_new.cXX.txt -ap -lmm -T 36
#mv -f Output/ vcf2gwas_oil_stderrs_stability

# phenotype stderr stability over all environments as input, no power transform
vcf2gwas -v ../GATK_SAM_population/IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf -pf oil_stderrs_stabilities.csv -k SAM_Kinship_new.cXX.txt -ap -lmm -T 36
mv -f Output/ vcf2gwas_oil_stderrs_stability_no_transform
