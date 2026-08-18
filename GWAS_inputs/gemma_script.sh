#!/bin/bash

# create kinship matrix
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p BC_oil_contents_GEMMA.csv -n 1 -gk 1 -o SAM_Kinship_new

# run multivariate GWAS on BC
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p BC_oil_contents_GEMMA.csv -k output/SAM_Kinship_new.cXX.txt -n 1 2 3 4 -lmm 1 -o GEMMA_BC

# run multivariate GWAS on GA
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p GA_oil_contents_GEMMA.csv -k output/SAM_Kinship_new.cXX.txt -n 1 2 3 4 -lmm 1 -o GEMMA_GA

# run multivariate GWAS on IA_2014
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p IA_2014_oil_contents_GEMMA.csv -k output/SAM_Kinship_new.cXX.txt -n 1 2 3 4 -lmm 1 -o GEMMA_IA_2014

# run multivariate GWAS on IA_2010
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p IA_oil_contents_GEMMA.csv -k output/SAM_Kinship_new.cXX.txt -n 1 2 3 4 -lmm 1 -o GEMMA_IA_2010

# run multivariate GWAS on IA_2013
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p IA_2013_oil_contents_GEMMA.csv -k output/SAM_Kinship_new.cXX.txt -n 1 2 3 4 -lmm 1 -o GEMMA_IA_2013

# run multivariate GWAS on MN_2015
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p MN_2015_oil_contents_GEMMA.csv -k output/SAM_Kinship_new.cXX.txt -n 1 2 3 4 -lmm 1 -o GEMMA_MN_2015

# run multivariate GWAS on MN_2016_early
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p MN_2016_day_1_oil_contents_GEMMA.csv -k output/SAM_Kinship_new.cXX.txt -n 1 2 3 4 -lmm 1 -o GEMMA_MN_2016_early

# run multivariate GWAS on MN_2016_late
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p MN_2016_day_2_oil_contents_GEMMA.csv -k output/SAM_Kinship_new.cXX.txt -n 1 2 3 4 -lmm 1 -o GEMMA_MN_2016_late

#########################################

# run multivariate GWAS on BC stderrs
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p BC_oil_contents_stderr_GEMMA.csv -k output/SAM_Kinship_new.cXX.txt -n 1 2 3 4 -lmm 1 -o GEMMA_BC_stderr

# run multivariate GWAS on GA stderrs
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p GA_oil_contents_stderr_GEMMA.csv -k output/SAM_Kinship_new.cXX.txt -n 1 2 3 4 -lmm 1 -o GEMMA_GA_stderr

# run multivariate GWAS on IA_2014 stderrs
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p IA_2014_oil_contents_stderr_GEMMA.csv -k output/SAM_Kinship_new.cXX.txt -n 1 2 3 4 -lmm 1 -o GEMMA_IA_2014_stderr

# run multivariate GWAS on IA_2013 stderrs
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p IA_2013_oil_contents_stderr_GEMMA.csv -k output/SAM_Kinship_new.cXX.txt -n 1 2 3 4 -lmm 1 -o GEMMA_IA_2013_stderr

# run multivariate GWAS on MN_2016_early stderrs
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p MN_2016_day_1_oil_contents_stderr_GEMMA.csv -k output/SAM_Kinship_new.cXX.txt -n 1 2 3 4 -lmm 1 -o GEMMA_MN_2016_early_stderr

# run multivariate GWAS on MN_2016_late stderrs
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p MN_2016_day_2_oil_contents_stderr_GEMMA.csv -k output/SAM_Kinship_new.cXX.txt -n 1 2 3 4 -lmm 1 -o GEMMA_MN_2016_late_stderr

########################################

# run multivariate GWAS on BC with FAD covariate
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p BC_oil_contents_GEMMA.csv -k output/SAM_Kinship_new.cXX.txt -n 1 2 3 4 -lmm 1 -c FAD2-1_mutant_status_GEMMA.csv -o GEMMA_BC_cov

# run multivariate GWAS on GA with FAD covariate
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p GA_oil_contents_GEMMA.csv -k output/SAM_Kinship_new.cXX.txt -n 1 2 3 4 -lmm 1 -c FAD2-1_mutant_status_GEMMA.csv -o GEMMA_GA_cov

# run multivariate GWAS on IA_2014 with FAD covariate
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p IA_2014_oil_contents_GEMMA.csv -k output/SAM_Kinship_new.cXX.txt -n 1 2 3 4 -lmm 1 -c FAD2-1_mutant_status_GEMMA.csv -o GEMMA_IA_2014_cov

# run multivariate GWAS on IA_2010 with FAD covariate
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p IA_oil_contents_GEMMA.csv -k output/SAM_Kinship_new.cXX.txt -n 1 2 3 4 -lmm 1 -c FAD2-1_mutant_status_GEMMA.csv -o GEMMA_IA_2010_cov

# run multivariate GWAS on IA_2013 with FAD covariate
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p IA_2013_oil_contents_GEMMA.csv -k output/SAM_Kinship_new.cXX.txt -n 1 2 3 4 -lmm 1 -c FAD2-1_mutant_status_GEMMA.csv -o GEMMA_IA_2013_cov

# run multivariate GWAS on MN_2015 with FAD covariate
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p MN_2015_oil_contents_GEMMA.csv -k output/SAM_Kinship_new.cXX.txt -n 1 2 3 4 -lmm 1 -c FAD2-1_mutant_status_GEMMA.csv -o GEMMA_MN_2015_cov

# run multivariate GWAS on MN_2016_early with FAD covariate
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p MN_2016_day_1_oil_contents_GEMMA.csv -k output/SAM_Kinship_new.cXX.txt -n 1 2 3 4 -lmm 1 -c FAD2-1_mutant_status_GEMMA.csv -o GEMMA_MN_2016_early_cov

# run multivariate GWAS on MN_2016_late with FAD covariate
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p MN_2016_day_2_oil_contents_GEMMA.csv -k output/SAM_Kinship_new.cXX.txt -n 1 2 3 4 -lmm 1 -c FAD2-1_mutant_status_GEMMA.csv -o GEMMA_MN_2016_late_cov

###########################################

# run multivariate GWAS on BC, omitting FAD mutants
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p BC_oil_contents_no_FAD_GEMMA.csv -k output/SAM_Kinship_new.cXX.txt -n 1 2 3 4 -lmm 1 -o GEMMA_BC_no_FAD

# run multivariate GWAS on GA, omitting FAD mutants
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p GA_oil_contents_no_FAD_GEMMA.csv -k output/SAM_Kinship_new.cXX.txt -n 1 2 3 4 -lmm 1 -o GEMMA_GA_no_FAD

# run multivariate GWAS on IA_2014, omitting FAD mutants
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p IA_2014_oil_contents_no_FAD_GEMMA.csv -k output/SAM_Kinship_new.cXX.txt -n 1 2 3 4 -lmm 1 -o GEMMA_IA_2014_no_FAD

# run multivariate GWAS on IA_2010, omitting FAD mutants
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p IA_oil_contents_no_FAD_GEMMA.csv -k output/SAM_Kinship_new.cXX.txt -n 1 2 3 4 -lmm 1 -o GEMMA_IA_2010_no_FAD

# run multivariate GWAS on IA_2013, omitting FAD mutants
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p IA_2013_oil_contents_no_FAD_GEMMA.csv -k output/SAM_Kinship_new.cXX.txt -n 1 2 3 4 -lmm 1 -o GEMMA_IA_2013_no_FAD

# run multivariate GWAS on MN_2015, omitting FAD mutants
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p MN_2015_oil_contents_no_FAD_GEMMA.csv -k output/SAM_Kinship_new.cXX.txt -n 1 2 3 4 -lmm 1 -o GEMMA_MN_2015_no_FAD

# run multivariate GWAS on MN_2016_early, omitting FAD mutants
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p MN_2016_day_1_oil_contents_no_FAD_GEMMA.csv -k output/SAM_Kinship_new.cXX.txt -n 1 2 3 4 -lmm 1 -o GEMMA_MN_2016_early_no_FAD

# run multivariate GWAS on MN_2016_late, omitting FAD mutants
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p MN_2016_day_2_oil_contents_no_FAD_GEMMA.csv -k output/SAM_Kinship_new.cXX.txt -n 1 2 3 4 -lmm 1 -o GEMMA_MN_2016_late_no_FAD

############################################

# run multivariate GWAS on stability across climates
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p oil_stabilities_GEMMA.csv -k output/SAM_Kinship_new.cXX.txt -n 1 2 3 4 -lmm 1 -o GEMMA_stability


# run multivariate GWAS on stability across climates with FAD covariate
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p oil_stabilities_GEMMA.csv -k output/SAM_Kinship_new.cXX.txt -n 1 2 3 4 -lmm 1 -c FAD2-1_mutant_status_GEMMA.csv -o GEMMA_stability_cov

# run multivariate GWAS on stability across climates without FAD samples
~/software/gemma-0.98.5-linux-static-AMD64 -bfile IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode -p oil_stabilities_no_FAD_GEMMA.csv -k output/SAM_Kinship_new.cXX.txt -n 1 2 3 4 -lmm 1 -o GEMMA_stability_no_FAD

exit 0