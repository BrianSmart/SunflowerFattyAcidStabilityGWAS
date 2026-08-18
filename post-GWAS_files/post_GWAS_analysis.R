# analysis of GEMMA output
library(qqman)
library(ggpubr)
library(dplyr)
#library(qqplotr)
library(stringr)
library(naniar)
library(readxl)
#library(arsenal)
library(CMplot)
library(tidyverse)
#library(xtable)
library(VennDiagram)
setwd("~/Documents/HDD/Sunflower/New_Plots")

# Horrible copy-pasting, look for a better way next time

#############################################
### GWAS on oil composition, no covariate ###
#############################################

# all environments, oleic
GWAS_BC_oleic <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_BC/Linear Mixed Model/Oleic/Oleic_20221110_130817/Oleic_mod_sub_BC_oil_contents.part3_IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.assoc.txt", header = T)
GWAS_GA_oleic <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_GA/Linear Mixed Model/Oleic/Oleic_20221110_133800/Oleic_mod_sub_GA_oil_contents.part3_IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.assoc.txt", header = T)
GWAS_IA_oleic <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_IA/Linear Mixed Model/Oleic/Oleic_20221110_141333/Oleic_mod_sub_IA_oil_contents.part3_IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.assoc.txt", header = T)
GWAS_MN_2015_oleic <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_MN_2015/Linear Mixed Model/Oleic/Oleic_20221110_143539/Oleic_mod_sub_MN_2015_oil_contents.part3_IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.assoc.txt", header = T)
GWAS_MN_2016_day_1_oleic <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_MN_2016_day_1/Linear Mixed Model/Oleic/Oleic_20221110_150538/Oleic_mod_sub_MN_2016_day_1_oil_contents.part3_IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.assoc.txt", header = T)
GWAS_MN_2016_day_2_oleic <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_MN_2016_day_2/Linear Mixed Model/Oleic/Oleic_20221110_153707/Oleic_mod_sub_MN_2016_day_2_oil_contents.part3_IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.assoc.txt", header = T)
GWAS_IA_2013_oleic <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_IA_2013/Linear Mixed Model/Oleic/Oleic_20221110_161713/Oleic_mod_sub_IA_2013_oil_contents.part3_IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.assoc.txt", header = T)
GWAS_IA_2014_oleic <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_IA_2014/Linear Mixed Model/Oleic/Oleic_20221110_165443/Oleic_mod_sub_IA_2014_oil_contents.part3_IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.assoc.txt", header = T)

# png("Manhattan_plot_vcf2gwas_all_envs_oleic.png", res=120, width = 1000, height = 4000)
# 
# layout(matrix(c(1,2,3,4,5,6,7,8), byrow = F),
#        heights =c(6,6,6,6,6,6,6,7))
# 
# manhattan(GWAS_BC_oleic, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("BC"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08), )
# manhattan(GWAS_GA_oleic, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("GA"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_oleic, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2015_oleic, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2015"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2016_day_1_oleic, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_day_1"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2016_day_2_oleic, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_day_2"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_2013_oleic, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA_2013"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_2014_oleic, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA_2014"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# dev.off()

# match the number of SNPs
GWAS_IA_oleic_temp <- merge(GWAS_BC_oleic,GWAS_IA_oleic, by = "rs", all.x = T, sort = F)

# circular manhattan plot
GWAS_univariate_oleic <- data.frame(SNP = GWAS_BC_oleic$rs, chr=GWAS_BC_oleic$chr, snp=GWAS_BC_oleic$ps, 
                            p_BC_oleic=GWAS_BC_oleic$p_wald, 
                            p_MN_2015_oleic=GWAS_MN_2015_oleic$p_wald,
                            p_MN_2016_early_oleic=GWAS_MN_2016_day_1_oleic$p_wald,
                            p_MN_2016_late_oleic=GWAS_MN_2016_day_2_oleic$p_wald,
                            p_IA_2010_oleic=GWAS_IA_oleic_temp$p_wald.y, 
                            p_IA_2013_oleic=GWAS_IA_2013_oleic$p_wald, 
                            p_IA_2014_oleic=GWAS_IA_2014_oleic$p_wald,
                            p_GA_oleic=GWAS_GA_oleic$p_wald)


CMplot(GWAS_univariate_oleic,type="p",plot.type="c",chr.labels=paste("Chr",c(1:17),sep=""),r=1,
       threshold=c(0.05/10937,2.26e-08), threshold.col=c("blue","red"),signal.line=1,signal.col=c("red","green"), 
       outward=F,cir.chr.h=0.5,chr.den.col="black",file="tiff", signal.cex = 0.6, amplify = T,
       file.output=T,verbose=TRUE,width=20,height=20)

# all environments, linoleic
GWAS_BC_linoleic <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_BC/Linear Mixed Model/Lino/Lino_20221110_130817/Lino_mod_sub_BC_oil_contents.part4_IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.assoc.txt", header = T)
GWAS_GA_linoleic <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_GA/Linear Mixed Model/Lino/Lino_20221110_133800/Lino_mod_sub_GA_oil_contents.part4_IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.assoc.txt", header = T)
GWAS_IA_linoleic <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_IA/Linear Mixed Model/Lino/Lino_20221110_141333/Lino_mod_sub_IA_oil_contents.part4_IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.assoc.txt", header = T)
GWAS_MN_2015_linoleic <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_MN_2015/Linear Mixed Model/Lino/Lino_20221110_143539/Lino_mod_sub_MN_2015_oil_contents.part4_IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.assoc.txt", header = T)
GWAS_MN_2016_day_1_linoleic <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_MN_2016_day_1/Linear Mixed Model/Lino/Lino_20221110_150538/Lino_mod_sub_MN_2016_day_1_oil_contents.part4_IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.assoc.txt", header = T)
GWAS_MN_2016_day_2_linoleic <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_MN_2016_day_2/Linear Mixed Model/Lino/Lino_20221110_153707/Lino_mod_sub_MN_2016_day_2_oil_contents.part4_IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.assoc.txt", header = T)
GWAS_IA_2013_linoleic <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_IA_2013/Linear Mixed Model/Lino/Lino_20221110_161713/Lino_mod_sub_IA_2013_oil_contents.part4_IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.assoc.txt", header = T)
GWAS_IA_2014_linoleic <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_IA_2014/Linear Mixed Model/Lino/Lino_20221110_165443/Lino_mod_sub_IA_2014_oil_contents.part4_IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.assoc.txt", header = T)
# 
# png("Manhattan_plot_vcf2gwas_all_envs_linoleic.png", res=120, width = 1000, height = 4000)
# 
# layout(matrix(c(1,2,3,4,5,6,7,8), byrow = F),
#        heights =c(6,6,6,6,6,6,6,7))
# 
# manhattan(GWAS_BC_linoleic, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("BC"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08), )
# manhattan(GWAS_GA_linoleic, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("GA"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_linoleic, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2015_linoleic, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2015"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2016_day_1_linoleic, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_day_1"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2016_day_2_linoleic, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_day_2"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_2013_linoleic, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA_2013"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_2014_linoleic, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA_2014"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# dev.off()

# match the number of SNPs
GWAS_IA_linoleic_temp <- merge(GWAS_BC_linoleic,GWAS_IA_linoleic, by = "rs", all.x = T, sort = F)

# circular manhattan plot
GWAS_univariate_linoleic <- data.frame(SNP = GWAS_BC_linoleic$rs, chr=GWAS_BC_linoleic$chr, snp=GWAS_BC_linoleic$ps, 
                                    p_BC_linoleic=GWAS_BC_linoleic$p_wald, 
                                    p_MN_2015_linoleic=GWAS_MN_2015_linoleic$p_wald,
                                    p_MN_2016_early_linoleic=GWAS_MN_2016_day_1_linoleic$p_wald,
                                    p_MN_2016_late_linoleic=GWAS_MN_2016_day_2_linoleic$p_wald,
                                    p_IA_2010_linoleic=GWAS_IA_linoleic_temp$p_wald.y, 
                                    p_IA_2013_linoleic=GWAS_IA_2013_linoleic$p_wald, 
                                    p_IA_2014_linoleic=GWAS_IA_2014_linoleic$p_wald,
                                    p_GA_linoleic=GWAS_GA_linoleic$p_wald)


CMplot(GWAS_univariate_linoleic,type="p",plot.type="c",chr.labels=paste("Chr",c(1:17),sep=""),r=1,
       threshold=c(0.05/10937,2.26e-08), threshold.col=c("blue","red"),signal.line=1,signal.col=c("red","green"), 
       outward=F,cir.chr.h=0.5,chr.den.col="black",file="tiff", signal.cex = 0.6, amplify = T,
       file.output=T,verbose=TRUE,width=20,height=20)

# all environments, palmitic
GWAS_BC_palm <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_BC/Linear Mixed Model/Palm/Palm_20221110_130817/Palm_mod_sub_BC_oil_contents.part1_IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.assoc.txt", header = T)
GWAS_GA_palm <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_GA/Linear Mixed Model/Palm/Palm_20221110_133800/Palm_mod_sub_GA_oil_contents.part1_IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.assoc.txt", header = T)
GWAS_IA_palm <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_IA/Linear Mixed Model/Palm/Palm_20221110_141333/Palm_mod_sub_IA_oil_contents.part1_IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.assoc.txt", header = T)
GWAS_MN_2015_palm <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_MN_2015/Linear Mixed Model/Palm/Palm_20221110_143539/Palm_mod_sub_MN_2015_oil_contents.part1_IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.assoc.txt", header = T)
GWAS_MN_2016_day_1_palm <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_MN_2016_day_1/Linear Mixed Model/Palm/Palm_20221110_150538/Palm_mod_sub_MN_2016_day_1_oil_contents.part1_IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.assoc.txt", header = T)
GWAS_MN_2016_day_2_palm <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_MN_2016_day_2/Linear Mixed Model/Palm/Palm_20221110_153707/Palm_mod_sub_MN_2016_day_2_oil_contents.part1_IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.assoc.txt", header = T)
GWAS_IA_2013_palm <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_IA_2013/Linear Mixed Model/Palm/Palm_20221110_161713/Palm_mod_sub_IA_2013_oil_contents.part1_IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.assoc.txt", header = T)
GWAS_IA_2014_palm <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_IA_2014/Linear Mixed Model/Palm/Palm_20221110_165443/Palm_mod_sub_IA_2014_oil_contents.part1_IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.assoc.txt", header = T)
# 
# png("Manhattan_plot_vcf2gwas_all_envs_palmitic.png", res=120, width = 1000, height = 4000)
# 
# layout(matrix(c(1,2,3,4,5,6,7,8), byrow = F),
#        heights =c(6,6,6,6,6,6,6,7))
# 
# manhattan(GWAS_BC_palm, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("BC"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08), )
# manhattan(GWAS_GA_palm, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("GA"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_palm, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2015_palm, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2015"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2016_day_1_palm, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_day_1"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2016_day_2_palm, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_day_2"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_2013_palm, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA_2013"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_2014_palm, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA_2014"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# dev.off()

# match the number of SNPs
GWAS_IA_palm_temp <- merge(GWAS_BC_palm,GWAS_IA_palm, by = "rs", all.x = T, sort = F)

# circular manhattan plot
GWAS_univariate_palm <- data.frame(SNP = GWAS_BC_palm$rs, chr=GWAS_BC_palm$chr, snp=GWAS_BC_palm$ps, 
                                    p_BC_palm=GWAS_BC_palm$p_wald, 
                                    p_MN_2015_palm=GWAS_MN_2015_palm$p_wald,
                                    p_MN_2016_early_palm=GWAS_MN_2016_day_1_palm$p_wald,
                                    p_MN_2016_late_palm=GWAS_MN_2016_day_2_palm$p_wald,
                                    p_IA_2010_palm=GWAS_IA_palm_temp$p_wald.y, 
                                    p_IA_2013_palm=GWAS_IA_2013_palm$p_wald, 
                                    p_IA_2014_palm=GWAS_IA_2014_palm$p_wald,
                                    p_GA_palm=GWAS_GA_palm$p_wald)


CMplot(GWAS_univariate_palm,type="p",plot.type="c",chr.labels=paste("Chr",c(1:17),sep=""),r=1,
       threshold=c(0.05/10937,2.26e-08), threshold.col=c("blue","red"),signal.line=1,signal.col=c("red","green"), 
       outward=F,cir.chr.h=0.5,chr.den.col="black",file="tiff", signal.cex = 0.6, amplify = T,
       file.output=T,verbose=TRUE,width=20,height=20)

# all environments, stearic
GWAS_BC_stea <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_BC/Linear Mixed Model/Stea/Stea_20221110_130817/Stea_mod_sub_BC_oil_contents.part2_IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.assoc.txt", header = T)
GWAS_GA_stea <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_GA/Linear Mixed Model/Stea/Stea_20221110_133800/Stea_mod_sub_GA_oil_contents.part2_IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.assoc.txt", header = T)
GWAS_IA_stea <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_IA/Linear Mixed Model/Stea/Stea_20221110_141333/Stea_mod_sub_IA_oil_contents.part2_IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.assoc.txt", header = T)
GWAS_MN_2015_stea <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_MN_2015/Linear Mixed Model/Stea/Stea_20221110_143539/Stea_mod_sub_MN_2015_oil_contents.part2_IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.assoc.txt", header = T)
GWAS_MN_2016_day_1_stea <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_MN_2016_day_1/Linear Mixed Model/Stea/Stea_20221110_150538/Stea_mod_sub_MN_2016_day_1_oil_contents.part2_IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.assoc.txt", header = T)
GWAS_MN_2016_day_2_stea <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_MN_2016_day_2/Linear Mixed Model/Stea/Stea_20221110_153707/Stea_mod_sub_MN_2016_day_2_oil_contents.part2_IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.assoc.txt", header = T)
GWAS_IA_2013_stea <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_IA_2013/Linear Mixed Model/Stea/Stea_20221110_161713/Stea_mod_sub_IA_2013_oil_contents.part2_IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.assoc.txt", header = T)
GWAS_IA_2014_stea <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_IA_2014/Linear Mixed Model/Stea/Stea_20221110_165443/Stea_mod_sub_IA_2014_oil_contents.part2_IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.assoc.txt", header = T)
# 
# png("Manhattan_plot_vcf2gwas_all_envs_stearic.png", res=120, width = 1000, height = 4000)
# 
# layout(matrix(c(1,2,3,4,5,6,7,8), byrow = F),
#        heights =c(6,6,6,6,6,6,6,7))
# 
# manhattan(GWAS_BC_stea, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("BC"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08), )
# manhattan(GWAS_GA_stea, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("GA"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_stea, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2015_stea, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2015"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2016_day_1_stea, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_day_1"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2016_day_2_stea, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_day_2"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_2013_stea, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA_2013"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_2014_stea, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA_2014"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# dev.off()

# match the number of SNPs
GWAS_IA_stea_temp <- merge(GWAS_BC_stea,GWAS_IA_stea, by = "rs", all.x = T, sort = F)

# circular manhattan plot
GWAS_univariate_stea <- data.frame(SNP = GWAS_BC_stea$rs, chr=GWAS_BC_stea$chr, snp=GWAS_BC_stea$ps, 
                                    p_BC_stea=GWAS_BC_stea$p_wald, 
                                    p_MN_2015_stea=GWAS_MN_2015_stea$p_wald,
                                    p_MN_2016_early_stea=GWAS_MN_2016_day_1_stea$p_wald,
                                    p_MN_2016_late_stea=GWAS_MN_2016_day_2_stea$p_wald,
                                    p_IA_2010_stea=GWAS_IA_stea_temp$p_wald.y, 
                                    p_IA_2013_stea=GWAS_IA_2013_stea$p_wald, 
                                    p_IA_2014_stea=GWAS_IA_2014_stea$p_wald,
                                    p_GA_stea=GWAS_GA_stea$p_wald)


CMplot(GWAS_univariate_stea,type="p",plot.type="c",chr.labels=paste("Chr",c(1:17),sep=""),r=1,
       threshold=c(0.05/10937,2.26e-08), threshold.col=c("blue","red"),signal.line=1,signal.col=c("red","green"), 
       outward=F,cir.chr.h=0.5,chr.den.col="black",file="tiff", signal.cex = 0.6, amplify = T,
       file.output=T,verbose=TRUE,width=20,height=20)

#################################################
### manhattan plots for FAD covariate in GWAS ###
#################################################
# all environments, oleic
GWAS_BC_oleic_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_BC_FAD_covariate/Linear Mixed Model/Oleic/Oleic_20221114_155122/assoc.txt", header = T)
GWAS_GA_oleic_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_GA_FAD_covariate/Linear Mixed Model/Oleic/Oleic_20221114_152735/assoc.txt", header = T)
GWAS_IA_oleic_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_IA_FAD_covariate/Linear Mixed Model/Oleic/Oleic_20221114_161432/assoc.txt", header = T)
GWAS_MN_2015_oleic_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_MN_2015_FAD_covariate/Linear Mixed Model/Oleic/Oleic_20221114_163513/assoc.txt", header = T)
GWAS_MN_2016_day_1_oleic_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_MN_2016_day_1_FAD_covariate/Linear Mixed Model/Oleic/Oleic_20221114_165847/assoc.txt", header = T)
GWAS_MN_2016_day_2_oleic_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_MN_2016_day_2_FAD_covariate/Linear Mixed Model/Oleic/Oleic_20221114_172513/assoc.txt", header = T)
GWAS_IA_2013_oleic_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_IA_2013_FAD_covariate/Linear Mixed Model/Oleic/Oleic_20221114_175316/assoc.txt", header = T)
GWAS_IA_2014_oleic_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_IA_2014_FAD_covariate/Linear Mixed Model/Oleic/Oleic_20221110_194946/assoc.txt", header = T)
# 
# png("Manhattan_plot_vcf2gwas_all_envs_oleic_FAD.png", res=120, width = 1000, height = 4000)
# 
# layout(matrix(c(1,2,3,4,5,6,7,8), byrow = F),
#        heights =c(6,6,6,6,6,6,6,7))
# 
# manhattan(GWAS_BC_oleic_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("BC"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08), )
# manhattan(GWAS_GA_oleic_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("GA"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_oleic_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2015_oleic_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2015"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2016_day_1_oleic_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_day_1"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2016_day_2_oleic_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_day_2"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_2013_oleic_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA_2013"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_2014_oleic_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA_2014"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# dev.off()

# match the number of SNPs
GWAS_IA_oleic_FAD_temp <- merge(GWAS_BC_oleic_FAD,GWAS_IA_oleic_FAD, by = "rs", all.x = T, sort = F)

# circular manhattan plot
GWAS_univariate_oleic_FAD <- data.frame(SNP = GWAS_BC_oleic_FAD$rs, chr=GWAS_BC_oleic_FAD$chr, snp=GWAS_BC_oleic_FAD$ps, 
                                    p_BC_oleic_FAD=GWAS_BC_oleic_FAD$p_wald, 
                                    p_MN_2015_oleic_FAD=GWAS_MN_2015_oleic_FAD$p_wald,
                                    p_MN_2016_early_oleic_FAD=GWAS_MN_2016_day_1_oleic_FAD$p_wald,
                                    p_MN_2016_late_oleic_FAD=GWAS_MN_2016_day_2_oleic_FAD$p_wald,
                                    p_IA_2010_oleic_FAD=GWAS_IA_oleic_FAD_temp$p_wald.y, 
                                    p_IA_2013_oleic_FAD=GWAS_IA_2013_oleic_FAD$p_wald, 
                                    p_IA_2014_oleic_FAD=GWAS_IA_2014_oleic_FAD$p_wald,
                                    p_GA_oleic_FAD=GWAS_GA_oleic_FAD$p_wald)


CMplot(GWAS_univariate_oleic_FAD,type="p",plot.type="c",chr.labels=paste("Chr",c(1:17),sep=""),r=1,
       threshold=c(0.05/10937,2.26e-08), threshold.col=c("blue","red"),signal.line=1,signal.col=c("red","green"), 
       outward=F,cir.chr.h=0.5,chr.den.col="black",file="tiff", signal.cex = 0.6, amplify = T,
       file.output=T,verbose=TRUE,width=20,height=20)

# all environments, linoleic
GWAS_BC_linoleic_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_BC_FAD_covariate/Linear Mixed Model/Lino/Lino_20221114_155122/assoc.txt", header = T)
GWAS_GA_linoleic_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_GA_FAD_covariate/Linear Mixed Model/Lino/Lino_20221114_152735/assoc.txt", header = T)
GWAS_IA_linoleic_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_IA_FAD_covariate/Linear Mixed Model/Lino/Lino_20221114_161432/assoc.txt", header = T)
GWAS_MN_2015_linoleic_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_MN_2015_FAD_covariate/Linear Mixed Model/Lino/Lino_20221114_163513/assoc.txt", header = T)
GWAS_MN_2016_day_1_linoleic_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_MN_2016_day_1_FAD_covariate/Linear Mixed Model/Lino/Lino_20221114_165847/assoc.txt", header = T)
GWAS_MN_2016_day_2_linoleic_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_MN_2016_day_2_FAD_covariate/Linear Mixed Model/Lino/Lino_20221114_172513/assoc.txt", header = T)
GWAS_IA_2013_linoleic_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_IA_2013_FAD_covariate/Linear Mixed Model/Lino/Lino_20221114_175316/assoc.txt", header = T)
GWAS_IA_2014_linoleic_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_IA_2014_FAD_covariate/Linear Mixed Model/Lino/Lino_20221110_194946/assoc.txt", header = T)
# 
# png("Manhattan_plot_vcf2gwas_all_envs_linoleic_FAD.png", res=120, width = 1000, height = 4000)
# 
# layout(matrix(c(1,2,3,4,5,6,7,8), byrow = F),
#        heights =c(6,6,6,6,6,6,6,7))
# 
# manhattan(GWAS_BC_linoleic_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("BC"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08), )
# manhattan(GWAS_GA_linoleic_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("GA"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_linoleic_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2015_linoleic_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2015"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2016_day_1_linoleic_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_day_1"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2016_day_2_linoleic_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_day_2"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_2013_linoleic_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA_2013"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_2014_linoleic_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA_2014"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# dev.off()

GWAS_IA_linoleic_FAD_temp <- merge(GWAS_BC_linoleic_FAD,GWAS_IA_linoleic_FAD, by = "rs", all.x = T, sort = F)

# circular manhattan plot
GWAS_univariate_linoleic_FAD <- data.frame(SNP = GWAS_BC_linoleic_FAD$rs, chr=GWAS_BC_linoleic_FAD$chr, snp=GWAS_BC_linoleic_FAD$ps, 
                                        p_BC_linoleic_FAD=GWAS_BC_linoleic_FAD$p_wald, 
                                        p_MN_2015_linoleic_FAD=GWAS_MN_2015_linoleic_FAD$p_wald,
                                        p_MN_2016_early_linoleic_FAD=GWAS_MN_2016_day_1_linoleic_FAD$p_wald,
                                        p_MN_2016_late_linoleic_FAD=GWAS_MN_2016_day_2_linoleic_FAD$p_wald,
                                        p_IA_2010_linoleic_FAD=GWAS_IA_linoleic_FAD_temp$p_wald.y, 
                                        p_IA_2013_linoleic_FAD=GWAS_IA_2013_linoleic_FAD$p_wald, 
                                        p_IA_2014_linoleic_FAD=GWAS_IA_2014_linoleic_FAD$p_wald,
                                        p_GA_linoleic_FAD=GWAS_GA_linoleic_FAD$p_wald)


CMplot(GWAS_univariate_linoleic_FAD,type="p",plot.type="c",chr.labels=paste("Chr",c(1:17),sep=""),r=1,
       threshold=c(0.05/10937,2.26e-08), threshold.col=c("blue","red"),signal.line=1,signal.col=c("red","green"), 
       outward=F,cir.chr.h=0.5,chr.den.col="black",file="tiff", signal.cex = 0.6, amplify = T,
       file.output=T,verbose=TRUE,width=20,height=20)

# all environments, palmitic
GWAS_BC_palm_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_BC_FAD_covariate/Linear Mixed Model/Palm/Palm_20221114_155122/assoc.txt", header = T)
GWAS_GA_palm_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_GA_FAD_covariate/Linear Mixed Model/Palm/Palm_20221114_152735/assoc.txt", header = T)
GWAS_IA_palm_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_IA_FAD_covariate/Linear Mixed Model/Palm/Palm_20221114_161432/assoc.txt", header = T)
GWAS_MN_2015_palm_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_MN_2015_FAD_covariate/Linear Mixed Model/Palm/Palm_20221114_163513/assoc.txt", header = T)
GWAS_MN_2016_day_1_palm_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_MN_2016_day_1_FAD_covariate/Linear Mixed Model/Palm/Palm_20221114_165847/assoc.txt", header = T)
GWAS_MN_2016_day_2_palm_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_MN_2016_day_2_FAD_covariate/Linear Mixed Model/Palm/Palm_20221114_172513/assoc.txt", header = T)
GWAS_IA_2013_palm_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_IA_2013_FAD_covariate/Linear Mixed Model/Palm/Palm_20221114_175316/assoc.txt", header = T)
GWAS_IA_2014_palm_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_IA_2014_FAD_covariate/Linear Mixed Model/Palm/Palm_20221110_194946/assoc.txt", header = T)
# 
# png("Manhattan_plot_/GEMMA_vcf2gwas_outputs/vcf2gwas_all_envs_palmitic_FAD.png", res=120, width = 1000, height = 4000)
# 
# layout(matrix(c(1,2,3,4,5,6,7,8), byrow = F),
#        heights =c(6,6,6,6,6,6,6,7))
# 
# manhattan(GWAS_BC_palm_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("BC"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08), )
# manhattan(GWAS_GA_palm_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("GA"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_palm_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2015_palm_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2015"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2016_day_1_palm_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_day_1"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2016_day_2_palm_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_day_2"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_2013_palm_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA_2013"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_2014_palm_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA_2014"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# dev.off()

GWAS_IA_palm_FAD_temp <- merge(GWAS_BC_palm_FAD,GWAS_IA_palm_FAD, by = "rs", all.x = T, sort = F)

# circular manhattan plot
GWAS_univariate_palm_FAD <- data.frame(SNP = GWAS_BC_palm_FAD$rs, chr=GWAS_BC_palm_FAD$chr, snp=GWAS_BC_palm_FAD$ps, 
                                        p_BC_palm_FAD=GWAS_BC_palm_FAD$p_wald, 
                                        p_MN_2015_palm_FAD=GWAS_MN_2015_palm_FAD$p_wald,
                                        p_MN_2016_early_palm_FAD=GWAS_MN_2016_day_1_palm_FAD$p_wald,
                                        p_MN_2016_late_palm_FAD=GWAS_MN_2016_day_2_palm_FAD$p_wald,
                                        p_IA_2010_palm_FAD=GWAS_IA_palm_FAD_temp$p_wald.y, 
                                        p_IA_2013_palm_FAD=GWAS_IA_2013_palm_FAD$p_wald, 
                                        p_IA_2014_palm_FAD=GWAS_IA_2014_palm_FAD$p_wald,
                                        p_GA_palm_FAD=GWAS_GA_palm_FAD$p_wald)


CMplot(GWAS_univariate_palm_FAD,type="p",plot.type="c",chr.labels=paste("Chr",c(1:17),sep=""),r=1,
       threshold=c(0.05/10937,2.26e-08), threshold.col=c("blue","red"),signal.line=1,signal.col=c("red","green"), 
       outward=F,cir.chr.h=0.5,chr.den.col="black",file="tiff", signal.cex = 0.6, amplify = T,
       file.output=T,verbose=TRUE,width=20,height=20)

# all environments, stearic
GWAS_BC_stea_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_BC_FAD_covariate/Linear Mixed Model/Stea/Stea_20221114_155122/assoc.txt", header = T)
GWAS_GA_stea_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_GA_FAD_covariate/Linear Mixed Model/Stea/Stea_20221114_152735/assoc.txt", header = T)
GWAS_IA_stea_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_IA_FAD_covariate/Linear Mixed Model/Stea/Stea_20221114_161432/assoc.txt", header = T)
GWAS_MN_2015_stea_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_MN_2015_FAD_covariate/Linear Mixed Model/Stea/Stea_20221114_163513/assoc.txt", header = T)
GWAS_MN_2016_day_1_stea_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_MN_2016_day_1_FAD_covariate/Linear Mixed Model/Stea/Stea_20221114_165847/assoc.txt", header = T)
GWAS_MN_2016_day_2_stea_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_MN_2016_day_2_FAD_covariate/Linear Mixed Model/Stea/Stea_20221114_172513/assoc.txt", header = T)
GWAS_IA_2013_stea_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_IA_2013_FAD_covariate/Linear Mixed Model/Stea/Stea_20221114_175316/assoc.txt", header = T)
GWAS_IA_2014_stea_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_IA_2014_FAD_covariate/Linear Mixed Model/Stea/Stea_20221110_194946/assoc.txt", header = T)
# 
# png("Manhattan_plot_/GEMMA_vcf2gwas_outputs/vcf2gwas_all_envs_stearic_FAD.png", res=120, width = 1000, height = 4000)
# 
# layout(matrix(c(1,2,3,4,5,6,7,8), byrow = F),
#        heights =c(6,6,6,6,6,6,6,7))
# 
# manhattan(GWAS_BC_stea_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("BC"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08), )
# manhattan(GWAS_GA_stea_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("GA"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_stea_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2015_stea_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2015"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2016_day_1_stea_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_day_1"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2016_day_2_stea_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_day_2"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_2013_stea_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA_2013"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_2014_stea_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA_2014"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# dev.off()

GWAS_IA_stea_FAD_temp <- merge(GWAS_BC_stea_FAD,GWAS_IA_stea_FAD, by = "rs", all.x = T, sort = F)

# circular manhattan plot
GWAS_univariate_stea_FAD <- data.frame(SNP = GWAS_BC_stea_FAD$rs, chr=GWAS_BC_stea_FAD$chr, snp=GWAS_BC_stea_FAD$ps, 
                                        p_BC_stea_FAD=GWAS_BC_stea_FAD$p_wald, 
                                        p_MN_2015_stea_FAD=GWAS_MN_2015_stea_FAD$p_wald,
                                        p_MN_2016_early_stea_FAD=GWAS_MN_2016_day_1_stea_FAD$p_wald,
                                        p_MN_2016_late_stea_FAD=GWAS_MN_2016_day_2_stea_FAD$p_wald,
                                        p_IA_2010_stea_FAD=GWAS_IA_stea_FAD_temp$p_wald.y, 
                                        p_IA_2013_stea_FAD=GWAS_IA_2013_stea_FAD$p_wald, 
                                        p_IA_2014_stea_FAD=GWAS_IA_2014_stea_FAD$p_wald,
                                        p_GA_stea_FAD=GWAS_GA_stea_FAD$p_wald)


CMplot(GWAS_univariate_stea_FAD,type="p",plot.type="c",chr.labels=paste("Chr",c(1:17),sep=""),r=1,
       threshold=c(0.05/10937,2.26e-08), threshold.col=c("blue","red"),signal.line=1,signal.col=c("red","green"), 
       outward=F,cir.chr.h=0.5,chr.den.col="black",file="tiff", signal.cex = 0.6, amplify = T,
       file.output=T,verbose=TRUE,width=20,height=20)

#######################################
### stability GWAS based on stderrs ###
#######################################

# all environments, oleic
GWAS_BC_oleic_stability <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_BC_stderr/Linear Mixed Model/Oleic_stderr/Oleic_stderr_20221110_173055/assoc.txt", header = T)
GWAS_GA_oleic_stability <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_GA_stderr/Linear Mixed Model/Oleic_stderr/Oleic_stderr_20221110_175105/assoc.txt", header = T)
GWAS_MN_2016_day_1_oleic_stability <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_MN_2016_day_1_stderr/Linear Mixed Model/Oleic_stderr/Oleic_stderr_20221110_181221/assoc.txt", header = T)
GWAS_MN_2016_day_2_oleic_stability <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_MN_2016_day_2_stderr/Linear Mixed Model/Oleic_stderr/Oleic_stderr_20221110_183406/assoc.txt", header = T)
GWAS_IA_2013_oleic_stability <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_IA_2013_stderr/Linear Mixed Model/Oleic_stderr/Oleic_stderr_20221110_190336/assoc.txt", header = T)
GWAS_IA_2014_oleic_stability <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_IA_2014_stderr/Linear Mixed Model/Oleic_stderr/Oleic_stderr_20221110_192518/assoc.txt", header = T)
#GWAS_oleic_stability_all <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_oil_stability/Linear Mixed Model/oleic_stability/oleic_stability_20221110_201935/assoc.txt", header = T)
#GWAS_oleic_stderrs <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_oil_stderrs_stability_no_transform/Linear Mixed Model/oleic_stderrs_stability/oleic_stderrs_stability_20221117_141215/assoc.txt")
# 
# png("Manhattan_plot_/GEMMA_vcf2gwas_outputs/vcf2gwas_all_envs_oleic_stability.png", res=120, width = 1000, height = 4000)
# 
# layout(matrix(c(1,2,3,4,5,6,7,8), byrow = F),
#        heights =c(6,6,6,6,6,6,6,7))
# 
# manhattan(GWAS_BC_oleic_stability, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("BC"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08), )
# manhattan(GWAS_GA_oleic_stability, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("GA"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2016_day_1_oleic_stability, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_day_1"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2016_day_2_oleic_stability, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_day_2"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_2013_oleic_stability, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA_2013"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_2014_oleic_stability, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA_2014"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_oleic_stability_all, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("Eberhart-Russell stability on oil content"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_oleic_stderrs, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("Eberhart-Russell stability on oil content stderr"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# 
# dev.off()

# match the number of SNPs
#GWAS_oleic_stability_all_temp <- merge(GWAS_BC_oleic_stability,GWAS_oleic_stability_all, by = "rs", all.x = T, sort = F)

# circular manhattan plot
GWAS_stability_oleic <- data.frame(SNP = GWAS_BC_oleic_stability$rs, chr=GWAS_BC_oleic_stability$chr, snp=GWAS_BC_oleic_stability$ps, 
                                    p_BC_oleic_stability=GWAS_BC_oleic_stability$p_wald, 
                                    p_MN_2016_early_oleic_stability=GWAS_MN_2016_day_1_oleic_stability$p_wald,
                                    p_MN_2016_late_oleic_stability=GWAS_MN_2016_day_2_oleic_stability$p_wald,
                                   p_IA_2013_oleic_stability=GWAS_IA_2013_oleic_stability$p_wald, 
                                   p_IA_2014_oleic_stability=GWAS_IA_2014_oleic_stability$p_wald,
                                   p_GA_oleic_stability=GWAS_GA_oleic_stability$p_wald 
)
                                   #p_all_envs=GWAS_oleic_stability_all_temp$p_wald.y)


CMplot(GWAS_stability_oleic,type="p",plot.type="c",chr.labels=paste("Chr",c(1:17),sep=""),r=1,
       threshold=c(0.05/10937,2.26e-08), threshold.col=c("blue","red"),signal.line=1,signal.col=c("red","green"), 
       outward=F,cir.chr.h=0.5,chr.den.col="black",file="tiff", signal.cex = 0.6, amplify = T,
       file.output=T,verbose=TRUE,width=20,height=20)

# all environments, linoleic
GWAS_BC_linoleic_stability <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_BC_stderr/Linear Mixed Model/Lino_stderr/Lino_stderr_20221110_173055/assoc.txt", header = T)
GWAS_GA_linoleic_stability <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_GA_stderr/Linear Mixed Model/Lino_stderr/Lino_stderr_20221110_175105/assoc.txt", header = T)
GWAS_MN_2016_day_1_linoleic_stability <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_MN_2016_day_1_stderr/Linear Mixed Model/Lino_stderr/Lino_stderr_20221110_181221/assoc.txt", header = T)
GWAS_MN_2016_day_2_linoleic_stability <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_MN_2016_day_2_stderr/Linear Mixed Model/Lino_stderr/Lino_stderr_20221110_183406/assoc.txt", header = T)
GWAS_IA_2013_linoleic_stability <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_IA_2013_stderr/Linear Mixed Model/Lino_stderr/Lino_stderr_20221110_190336/assoc.txt", header = T)
GWAS_IA_2014_linoleic_stability <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_IA_2014_stderr/Linear Mixed Model/Lino_stderr/Lino_stderr_20221110_192518/assoc.txt", header = T)
#GWAS_linoleic_stability_all <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_oil_stability/Linear Mixed Model/linoleic_stability/linoleic_stability_20221110_201935/assoc.txt", header = T)
#GWAS_linoleic_stderrs <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_oil_stderrs_stability_no_transform/Linear Mixed Model/linoleic_stderrs_stability/linoleic_stderrs_stability_20221117_141215/assoc.txt")
# 
# png("Manhattan_plot_/GEMMA_vcf2gwas_outputs/vcf2gwas_all_envs_linoleic_stability.png", res=120, width = 1000, height = 4000)
# 
# layout(matrix(c(1,2,3,4,5,6,7,8), byrow = F),
#        heights =c(6,6,6,6,6,6,6,7))
# 
# manhattan(GWAS_BC_linoleic_stability, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("BC"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08), )
# manhattan(GWAS_GA_linoleic_stability, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("GA"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2016_day_1_linoleic_stability, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_day_1"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2016_day_2_linoleic_stability, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_day_2"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_2013_linoleic_stability, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA_2013"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_2014_linoleic_stability, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA_2014"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_linoleic_stability_all, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("Eberhart-Russell stability on oil content"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_linoleic_stderrs, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("Eberhart-Russell stability on oil content stderr"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# 
# dev.off()

GWAS_stability_linoleic <- data.frame(SNP = GWAS_BC_linoleic_stability$rs, chr=GWAS_BC_linoleic_stability$chr, snp=GWAS_BC_linoleic_stability$ps, 
                                   p_BC_linoleic_stability=GWAS_BC_linoleic_stability$p_wald, 
                                   p_MN_2016_early_linoleic_stability=GWAS_MN_2016_day_1_linoleic_stability$p_wald,
                                   p_MN_2016_late_linoleic_stability=GWAS_MN_2016_day_2_linoleic_stability$p_wald,
                                   p_IA_2013_linoleic_stability=GWAS_IA_2013_linoleic_stability$p_wald, 
                                   p_IA_2014_linoleic_stability=GWAS_IA_2014_linoleic_stability$p_wald,
                                   p_GA_linoleic_stability=GWAS_GA_linoleic_stability$p_wald
)
#p_all_envs=GWAS_linoleic_stability_all_temp$p_wald.y)


CMplot(GWAS_stability_linoleic,type="p",plot.type="c",chr.labels=paste("Chr",c(1:17),sep=""),r=1,
       threshold=c(0.05/10937,2.26e-08), threshold.col=c("blue","red"),signal.line=1,signal.col=c("red","green"), 
       outward=F,cir.chr.h=0.5,chr.den.col="black",file="tiff", signal.cex = 0.6, amplify = T,
       file.output=T,verbose=TRUE,width=20,height=20)

# all environments, palmitic
GWAS_BC_palm_stability <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_BC_stderr/Linear Mixed Model/Palm_stderr/Palm_stderr_20221110_173055/assoc.txt", header = T)
GWAS_GA_palm_stability <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_GA_stderr/Linear Mixed Model/Palm_stderr/Palm_stderr_20221110_175105/assoc.txt", header = T)
GWAS_MN_2016_day_1_palm_stability <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_MN_2016_day_1_stderr/Linear Mixed Model/Palm_stderr/Palm_stderr_20221110_181221/assoc.txt", header = T)
GWAS_MN_2016_day_2_palm_stability <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_MN_2016_day_2_stderr/Linear Mixed Model/Palm_stderr/Palm_stderr_20221110_183406/assoc.txt", header = T)
GWAS_IA_2013_palm_stability <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_IA_2013_stderr/Linear Mixed Model/Palm_stderr/Palm_stderr_20221110_190336/assoc.txt", header = T)
GWAS_IA_2014_palm_stability <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_IA_2014_stderr/Linear Mixed Model/Palm_stderr/Palm_stderr_20221110_192518/assoc.txt", header = T)
#GWAS_palm_stability_all <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_oil_stability/Linear Mixed Model/palmitic_stability/palmitic_stability_20221110_201935/assoc.txt", header = T)
#GWAS_palm_stderrs <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_oil_stderrs_stability_no_transform/Linear Mixed Model/palmitic_stderrs_stability/palmitic_stderrs_stability_20221117_141215/assoc.txt")
# 
# png("Manhattan_plot_/GEMMA_vcf2gwas_outputs/vcf2gwas_all_envs_palmitic_stability.png", res=120, width = 1000, height = 4000)
# 
# layout(matrix(c(1,2,3,4,5,6,7,8), byrow = F),
#        heights =c(6,6,6,6,6,6,6,7))
# 
# manhattan(GWAS_BC_palm_stability, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("BC"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08), )
# manhattan(GWAS_GA_palm_stability, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("GA"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2016_day_1_palm_stability, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_day_1"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2016_day_2_palm_stability, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_day_2"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_2013_palm_stability, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA_2013"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_2014_palm_stability, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA_2014"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_palm_stability_all, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("Eberhart-Russell stability on oil content"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_palm_stderrs, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("Eberhart-Russell stability on oil content stderr"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# 
# dev.off()

GWAS_stability_palm <- data.frame(SNP = GWAS_BC_palm_stability$rs, chr=GWAS_BC_palm_stability$chr, snp=GWAS_BC_palm_stability$ps, 
                                   p_BC_palm_stability=GWAS_BC_palm_stability$p_wald, 
                                   p_MN_2016_early_palm_stability=GWAS_MN_2016_day_1_palm_stability$p_wald,
                                   p_MN_2016_late_palm_stability=GWAS_MN_2016_day_2_palm_stability$p_wald,
                                   p_IA_2013_palm_stability=GWAS_IA_2013_palm_stability$p_wald, 
                                   p_IA_2014_palm_stability=GWAS_IA_2014_palm_stability$p_wald,
                                   p_GA_palm_stability=GWAS_GA_palm_stability$p_wald
)
#p_all_envs=GWAS_palm_stability_all_temp$p_wald.y)

CMplot(GWAS_stability_palm,type="p",plot.type="c",chr.labels=paste("Chr",c(1:17),sep=""),r=1,
       threshold=c(0.05/10937,2.26e-08), threshold.col=c("blue","red"),signal.line=1,signal.col=c("red","green"), 
       outward=F,cir.chr.h=0.5,chr.den.col="black",file="tiff", signal.cex = 0.6, amplify = T,
       file.output=T,verbose=TRUE,width=20,height=20)

# all environments, stearic
GWAS_BC_stea_stability <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_BC_stderr/Linear Mixed Model/Stea_stderr/Stea_stderr_20221110_173055/assoc.txt", header = T)
GWAS_GA_stea_stability <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_GA_stderr/Linear Mixed Model/Stea_stderr/Stea_stderr_20221110_175105/assoc.txt", header = T)
GWAS_MN_2016_day_1_stea_stability <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_MN_2016_day_1_stderr/Linear Mixed Model/Stea_stderr/Stea_stderr_20221110_181221/assoc.txt", header = T)
GWAS_MN_2016_day_2_stea_stability <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_MN_2016_day_2_stderr/Linear Mixed Model/Stea_stderr/Stea_stderr_20221110_183406/assoc.txt", header = T)
GWAS_IA_2013_stea_stability <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_IA_2013_stderr/Linear Mixed Model/Stea_stderr/Stea_stderr_20221110_190336/assoc.txt", header = T)
GWAS_IA_2014_stea_stability <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_IA_2014_stderr/Linear Mixed Model/Stea_stderr/Stea_stderr_20221110_192518/assoc.txt", header = T)
#GWAS_stea_stability_all <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_oil_stability/Linear Mixed Model/stearic_stability/stearic_stability_20221110_201935/assoc.txt", header = T)
#GWAS_stea_stderrs <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_oil_stderrs_stability_no_transform/Linear Mixed Model/stearic_stderrs_stability/stearic_stderrs_stability_20221117_141215/assoc.txt")
# 
# png("Manhattan_plot_/GEMMA_vcf2gwas_outputs/vcf2gwas_all_envs_stearic_stability.png", res=120, width = 1000, height = 4000)
# 
# layout(matrix(c(1,2,3,4,5,6,7,8), byrow = F),
#        heights =c(6,6,6,6,6,6,6,7))
# 
# manhattan(GWAS_BC_stea_stability, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("BC"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08), )
# manhattan(GWAS_GA_stea_stability, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("GA"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2016_day_1_stea_stability, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_day_1"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2016_day_2_stea_stability, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_day_2"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_2013_stea_stability, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA_2013"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_2014_stea_stability, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA_2014"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_stea_stability_all, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("Eberhart-Russell stability on oil content"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_stea_stderrs, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("Eberhart-Russell stability on oil content stderr"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# 
# dev.off()

GWAS_stability_stea <- data.frame(SNP = GWAS_BC_stea_stability$rs, chr=GWAS_BC_stea_stability$chr, snp=GWAS_BC_stea_stability$ps, 
                                   p_BC_stea_stability=GWAS_BC_stea_stability$p_wald, 
                                   p_MN_2016_early_stea_stability=GWAS_MN_2016_day_1_stea_stability$p_wald,
                                   p_MN_2016_late_stea_stability=GWAS_MN_2016_day_2_stea_stability$p_wald,
                                   p_IA_2013_stea_stability=GWAS_IA_2013_stea_stability$p_wald, 
                                   p_IA_2014_stea_stability=GWAS_IA_2014_stea_stability$p_wald,
                                   p_GA_stea_stability=GWAS_GA_stea_stability$p_wald
)
#p_all_envs=GWAS_stea_stability_all_temp$p_wald.y)


CMplot(GWAS_stability_stea,type="p",plot.type="c",chr.labels=paste("Chr",c(1:17),sep=""),r=1,
       threshold=c(0.05/10937,2.26e-08), threshold.col=c("blue","red"),signal.line=1,signal.col=c("red","green"), 
       outward=F,cir.chr.h=0.5,chr.den.col="black",file="tiff", signal.cex = 0.6, amplify = T,
       file.output=T,verbose=TRUE,width=20,height=20)

##################################################
### IA_2014 multivariate/univariate comparison ###
##################################################

GEMMA_IA_2014_oleic <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_univariate_oleic_IA_2014.assoc.txt", header = T)
GEMMA_IA_2014_oleic_no_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/vcf2gwas_IA_2014_no_FAD_mutants/Linear Mixed Model/Oleic/Oleic_20221114_125906/Oleic_mod_sub_IA_2014_oil_contents_no_FAD.part3_IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.assoc.txt")
GEMMA_IA_2014_palm <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_univariate_palm_IA_2014.assoc.txt", header = T)
GEMMA_IA_2014_stea <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_univariate_stea_IA_2014.assoc.txt", header = T)
GEMMA_IA_2014_lino <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_univariate_lino_IA_2014.assoc.txt", header = T)
GEMMA_IA_2014_multivariate <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_multivariate_IA_2014.assoc.txt", header = T)
#GEMMA_IA_2014_no_lino <- read.delim("../GEMMA_IA_2014_no_linoleic.assoc.txt", header = T)
GEMMA_IA_2014_no_FAD_multivariate <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_IA_2014_no_FAD.assoc.txt", header = T)
#GEMMA_IA_2014_FAD_covariate_multivariate <- read.delim("../GEMMA_vcf2gwas_outputs/IA_2014_FAD_cov.assoc.txt", header = T)
GEMMA_IA_2014_multivariate_2_PCs <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_IA_2014_2_PCs.assoc.txt", header = T)
GEMMA_IA_2014_no_FAD_multivariate_2_PCs <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_IA_2014_no_FAD_2_PCs.assoc.txt", header = T)
# 
# png("Manhattan_plot_GEMMA_single+multivariate_IA2014.png", res=100, width = 1600, height = 5500)
# layout(matrix(c(1,2,3,4,5,6,7), byrow = F),
#        heights =c(6,6,6,6,6,6,7))
# 
# manhattan(GEMMA_IA_2014_oleic, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("oleic"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08), )
# manhattan(GEMMA_IA_2014_lino, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("linoleic"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08), )
# manhattan(GEMMA_IA_2014_palm, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("palmitic"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08), )
# manhattan(GEMMA_IA_2014_stea, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("stearic"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08), )
# manhattan(GEMMA_IA_2014_multivariate, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("multivariate"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08), )
# #manhattan(GEMMA_IA_2014_no_lino, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("multivariate without linoleic"),
# #          suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08), )
# manhattan(GEMMA_IA_2014_no_FAD_multivariate, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("multivariate no FAD mutants"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08), )
# manhattan(GEMMA_IA_2014_FAD_covariate_multivariate, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("multivariate with FAD covariate"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08), )
# 
# dev.off()
# 
# # IA_2014 multivariate with and without FAD, 2 PC covariate
# png("Manhattan_plot_multivariate_IA_2014+2_PCs.png", res=100, width = 1600, height = 2000)
# layout(matrix(c(1,2,3,4), byrow = F),
#        heights =c(6,6,6,7))
# 
# manhattan(GEMMA_IA_2014_multivariate, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("multivariate all, no PCs"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08), )
# manhattan(GEMMA_IA_2014_multivariate_2_PCs, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("multivariate all, 2 PCs"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08), )
# manhattan(GEMMA_IA_2014_no_FAD_multivariate, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("multivariate no FAD, no PCs"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08), )
# manhattan(GEMMA_IA_2014_no_FAD_multivariate_2_PCs, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("multivariate no FAD, 2 PCs"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08), )
# dev.off()
# 
# ### create a qq plot to compare all sample GWAS with FAD mutant omission and univariate ###
# 
# # match df lengths
# GEMMA_IA_2014_no_FAD_multivariate_temp <- merge(GEMMA_IA_2014_multivariate, GEMMA_IA_2014_no_FAD_multivariate, by = c("chr", "ps"), all.x = T, sort = F)
# GEMMA_IA_2014_no_FAD_oleic_temp <- merge(GEMMA_IA_2014_multivariate, GEMMA_IA_2014_oleic_no_FAD, by = c("chr", "ps"), all.x = T, sort = F)
# 
# 
# GWAS_qqplots <- data.frame(name = rownames(GEMMA_IA_2014_multivariate), chr=GEMMA_IA_2014_multivariate$chr, snp=GEMMA_IA_2014_multivariate$ps, 
#                            p_IA_2014_multivariate=GEMMA_IA_2014_multivariate$p_wald, 
#                            p_IA_2014_univariate_oleic=GEMMA_IA_2014_oleic$p_wald,
#                            p_IA_2014_multivariate_no_FAD=GEMMA_IA_2014_no_FAD_multivariate_temp$p_wald.y,
#                            p_IA_2014_univariate_oleic_no_FAD=GEMMA_IA_2014_no_FAD_oleic_temp$p_wald.y)
# 
# CMplot(GWAS_qqplots,plot.type="q",col=c("dodgerblue1", "darkgoldenrod1", "olivedrab3", "brown2"),threshold=2.26e-08,
#        ylab.pos=2,signal.pch=c(6,4),signal.cex=1.2,signal.col="red",conf.int=TRUE,box=FALSE,multracks=
#          TRUE,cex.axis=2,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,width=10,height=10)
# 
# 
# 
# # manhattan of multivariate IA_2014
# png("Manhattan_plot_GEMMA_multivariate_IA_2014.png", res=100, width = 1600, height = 500)
# manhattan(GEMMA_IA_2014_multivariate, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", #main = ("multivariate"), 
#           suggestiveline = F, 
#           genomewideline = -log10(2.26e-08))
# dev.off()
# 
# # try to use a circos plot instead
# # create input dataframe
# IA_2014_input <- data.frame(SNP = GEMMA_IA_2014_oleic$rs, chr=GEMMA_IA_2014_oleic$chr, snp=GEMMA_IA_2014_oleic$ps, 
#                             p_stea=GEMMA_IA_2014_stea$p_wald, 
#                             p_palm=GEMMA_IA_2014_palm$p_wald,
#                             p_lino=GEMMA_IA_2014_lino$p_wald, 
#                             p_oleic=GEMMA_IA_2014_oleic$p_wald, 
#                             p_multi=GEMMA_IA_2014_multivariate$p_wald)
# 
# CMplot(IA_2014_input,type="p",plot.type="c",chr.labels=paste("Chr",c(1:17),sep=""),r=0.4,cir.legend=TRUE,
#        threshold=c(2.26e-08), threshold.col=c("blue"),signal.line=1,signal.col=c("red","green"), 
#        outward=F,cir.legend.col="black",cir.chr.h=0.5,chr.den.col="black",file="jpg", signal.cex = 0.1,
#        memo="",dpi=300,file.output=T,verbose=TRUE,width=20,height=20)
# 
# 
# # plot p-value differences between the two methods of dealing with FAD mutants
# GEMMA_IA_2014_FAD_multivariate_merge <- merge(GEMMA_IA_2014_no_FAD_multivariate, GEMMA_IA_2014_FAD_covariate_multivariate, by=c("chr", "ps"))
# GEMMA_IA_2014_FAD_multivariate_merge$p_val_diff <- abs(-log10(GEMMA_IA_2014_FAD_multivariate_merge$p_wald.x) - (-log10(GEMMA_IA_2014_FAD_multivariate_merge$p_wald.y)))
# 
# png("Manhattan_plot_multivariate_FAD_differences.png", res=100, width = 1600, height = 1200)
# manhattan(GEMMA_IA_2014_FAD_multivariate_merge, chr = "chr", bp = "ps", snp = "rs.x", p = "p_val_diff", main = ("Difference between FAD methods"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08), logp = F, ylab = "p_val_diff")
# dev.off()
# 
# 
# # look for the SNPs around the fad2-1 gene
# # fad2 = HanXRQr2_Chr14g0657451, pos 166019672 on chr 14 according to HAN412_Eugene_curated_v1_1.gff3
# GWAS_chr14 <- subset(GEMMA_IA_2014_oleic, chr == 14)
# gene_pos <- 166019672
# SNPs_around_fad2 <- subset(GWAS_chr14, ps > gene_pos - 1000000 & ps < gene_pos + 1000000)
# 
# 
# png("Manhattan_plot_IA_2014_fad2_marked.png", res=100, width = 1600, height = 1200)
# manhattan(subset(GEMMA_IA_2014_oleic, chr == 14), chr = "chr", bp = "ps", snp = "rs", p = "p_wald",
#           suggestiveline = -log10(0.00001), genomewideline = -log10(2.26e-08), xlim = c(100000000, 180000000),
#           highlight = SNPs_around_fad2$rs)
# dev.off()
# 

# # calculate lambda to quantify p-value inflation 
# sum(is.na(unlist(GWAS$p_wald))) # 0
# 
# lambda <- round(median(qchisq(GWAS$p_wald,1, lower.tail=FALSE)) / qchisq(0.5, 1), 3)
# 
# # qq plot
#png("qqplot_BC_vcf2gwas.png", res=100, width = 1000, height = 1000)
#png("qqplot_BC_GEMMA.png", res=100, width = 1000, height = 1000)
# png("qqplot_MN_2016_day_1_FAD_covariate_vcf2gwas.png", res=100, width = 1000, height = 1000)
# qq(GWAS$p_wald)
# #text(x= 1, y = 8, labels = paste0("lambda = ",lambda))
# 
# dev.off()
# 
# GWAS_MN_2016_day_1_no_FAD <- read.delim("../vcf2gwas_MN_2016_day_1_no_FAD_mutants/Linear Mixed Model/Palm+Stea+Oleic+Lino/Palm+Stea+Oleic+Lino_20221006_142958/Palm+Stea+Oleic+Lino_mod_sub_MN_2016_day_1_oil_contents_no_FAD.txt", header = T)
# 
# png("Manhattan_plot_MN_2016_day_1_no_FAD_mutants_vcf2gwas.png", res=100, width = 1600, height = 1200)
# manhattan_MN_2016_day_1 <- manhattan(GWAS_MN_2016_day_1_no_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_day_1 no FAD mutants"),
#                                      suggestiveline = -log10(0.00001), genomewideline = -log10(2.26e-08))
# dev.off()
# 
# 
# png("Manhattan_plot_MN_2016_day_1_no_FAD_mutants_vcf2gwas.png", res=120, width = 1600, height = 1400)
# layout(matrix(c(1,2), byrow = F),
#        heights =c(6,7))
# 
# manhattan(GWAS_IA_2013_stability, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("GWAS_IA_2013_stability"),
#           suggestiveline = -log10(0.00001), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_2014_stability, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("GWAS_IA_2014_stability"),
#           suggestiveline = -log10(0.00001), genomewideline = -log10(2.26e-08))
# dev.off()
# 
# 
# # format GEMMA result for input in LDBlockShow
# LDBlockShow_input <- GWAS_IA_2014[c(1,3,12)]
# write.table(LDBlockShow_input, file = "LDBlockShow_input.txt", sep = " ", row.names = F, col.names = F)
# 

#########################
### multivariate GWAS ###
#########################

# beta stability

GWAS_stability_multivariate <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_stability.assoc.txt", header = T)
#GWAS_stability_multivariate_norm <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_stability_norm.assoc.txt", header = T)
#GWAS_stability_multivariate_cov <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_stability_cov.assoc.txt", header = T)
GWAS_stability_multivariate_no_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_stability_no_FAD.assoc.txt", header = T)

png("Manhattan_plot_GEMMA_multivariate_stability.png", res=100, width = 1600, height = 1000)

layout(matrix(c(1,2), byrow = F),
       heights =c(6,7))

manhattan(GWAS_stability_multivariate, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("stability multivariate"),
          suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08), )
# manhattan(GWAS_stability_multivariate_norm, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("stability multivariate normalized"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# dev.off()
# 
# manhattan(GWAS_stability_multivariate_cov, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("stability multivariate + FAD covariate"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# 
manhattan(GWAS_stability_multivariate_no_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("stability multivariate, FAD mutants omitted"),
          suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
dev.off()

tiff("Manhattan_plot_GEMMA_multivariate_stability.tiff", width = 1500, height = 500)

manhattan(GWAS_stability_multivariate, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", #main = ("stability multivariate"),
          suggestiveline = F, genomewideline = -log10(2.26e-08), )
dev.off()

# mean FA contents
GWAS_BC_multivariate <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_BC.assoc.txt", header = T)
GWAS_GA_multivariate <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_GA.assoc.txt", header = T)
GWAS_IA_multivariate <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_IA_2010.assoc.txt", header = T)
GWAS_MN_2015_multivariate <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_MN_2015.assoc.txt", header = T)
GWAS_MN_2016_day_1_multivariate <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_MN_2016_early.assoc.txt", header = T)
GWAS_MN_2016_day_2_multivariate <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_MN_2016_late.assoc.txt", header = T)
GWAS_IA_2013_multivariate <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_IA_2013.assoc.txt", header = T)
GWAS_IA_2014_multivariate <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_IA_2014.assoc.txt", header = T)
# 
# png("Manhattan_plot_GEMMA_all_envs_multivariate.png", res=120, width = 1500, height = 4000)
# 
# layout(matrix(c(1,2,3,4,5,6,7,8), byrow = F),
#        heights =c(6,6,6,6,6,6,6,7))
# 
# manhattan(GWAS_BC_multivariate, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("BC"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08), )
# manhattan(GWAS_GA_multivariate, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("GA"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_multivariate, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2015_multivariate, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2015"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2016_day_1_multivariate, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_day_1"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2016_day_2_multivariate, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_day_2"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_2013_multivariate, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA_2013"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_2014_multivariate, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA_2014"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# dev.off()

# calculate FDR threshold. Too small for the multivariate p-value inflation!
#CalcThreshold(GWAS_BC_multivariate, sig.level = 0.05, method = "fdr")

# match df lengths
GWAS_IA_multivariate_temp <- merge(GWAS_BC_multivariate,GWAS_IA_multivariate, by = c("chr", "ps"), all.x = T, sort = F)
GWAS_stability_multivariate_temp <- merge(GWAS_BC_multivariate,GWAS_stability_multivariate, by = c("chr", "ps"), all.x = T, sort = F)

# create input df for CMplot
GWAS_multivariate <- data.frame(name = rownames(GWAS_BC_multivariate), chr=GWAS_BC_multivariate$chr, snp=GWAS_BC_multivariate$ps, 
                                p_BC=GWAS_BC_multivariate$p_wald, 
                                p_MN_2015=GWAS_MN_2015_multivariate$p_wald,
                                p_MN_2016_early=GWAS_MN_2016_day_1_multivariate$p_wald, 
                                p_MN_2016_late=GWAS_MN_2016_day_2_multivariate$p_wald,
                                p_IA_2010=GWAS_IA_multivariate_temp$p_wald.y, 
                                p_IA_2013=GWAS_IA_2013_multivariate$p_wald, 
                                p_IA_2014=GWAS_IA_2014_multivariate$p_wald,
                                p_GA=GWAS_GA_multivariate$p_wald,
                                p_stability=GWAS_stability_multivariate_temp$p_wald.y
                                )

# circular plot of multivariate GWAS
CMplot(GWAS_multivariate,type="p",plot.type="c",chr.labels=paste("Chr",c(1:17),sep=""),r=1,
       threshold=c(0.05/10937,2.26e-08), threshold.col=c("blue","red"),signal.line=1,signal.col=c("red","green"), 
       outward=F,cir.chr.h=0.5,chr.den.col="black",file="tiff", signal.cex = 0.6, amplify = T,
       file.output=T,verbose=TRUE,width=20,height=20)


# #FAD covariate
# # something is wrong with GWAS_MN_2016_day_1_multivariate_cov!!!
# GWAS_BC_multivariate_cov <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_BC_cov.assoc.txt", header = T)
# GWAS_GA_multivariate_cov <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_GA_cov.assoc.txt", header = T)
# GWAS_IA_multivariate_cov <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_IA_2010_cov.assoc.txt", header = T)
# GWAS_MN_2015_multivariate_cov <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_MN_2015_cov.assoc.txt", header = T)
# GWAS_MN_2016_day_1_multivariate_cov <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_MN_2016_early_cov.assoc.txt", header = T)
# GWAS_MN_2016_day_2_multivariate_cov <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_MN_2016_late_cov.assoc.txt", header = T)
# GWAS_IA_2013_multivariate_cov <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_IA_2013_cov.assoc.txt", header = T)
# GWAS_IA_2014_multivariate_cov <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_IA_2014_cov.assoc.txt", header = T)
# 
# png("Manhattan_plot_GEMMA_all_envs_multivariate_cov.png", res=120, width = 1500, height = 4000)
# 
# layout(matrix(c(1,2,3,4,5,6,7,8), byrow = F),
#        heights =c(6,6,6,6,6,6,6,7))
# 
# manhattan(GWAS_BC_multivariate_cov, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("BC"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08), )
# manhattan(GWAS_GA_multivariate_cov, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("GA"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_multivariate_cov, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2015_multivariate_cov, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2015"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2016_day_1_multivariate_cov, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_day_1"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2016_day_2_multivariate_cov, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_day_2"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_2013_multivariate_cov, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA_2013"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_2014_multivariate_cov, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA_2014"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# dev.off()
# 
# # match df lengths
# GWAS_IA_multivariate_cov_temp <- merge(GWAS_BC_multivariate_cov,GWAS_IA_multivariate_cov, by = c("chr", "ps"), all.x = T, sort = F)
# GWAS_stability_multivariate_cov_temp <- merge(GWAS_BC_multivariate_cov,GWAS_stability_multivariate_cov, by = c("chr", "ps"), all.x = T, sort = F)
# 
# 
# GWAS_multivariate_cov <- data.frame(name = rownames(GWAS_BC_multivariate_cov), chr=GWAS_BC_multivariate_cov$chr, snp=GWAS_BC_multivariate_cov$ps, 
#                                 p_cov_BC=GWAS_BC_multivariate_cov$p_wald, 
#                                 p_cov_MN_2015=GWAS_MN_2015_multivariate_cov$p_wald,
#                                 #p_MN_2016_early=GWAS_MN_2016_day_1_multivariate_cov$p_wald, # omitted for extremely high p-values 
#                                 p_cov_MN_2016_late=GWAS_MN_2016_day_2_multivariate_cov$p_wald, 
#                                 p_cov_IA_2010=GWAS_IA_multivariate_cov_temp$p_wald.y, 
#                                 p_cov_IA_2013=GWAS_IA_2013_multivariate_cov$p_wald, 
#                                 p_cov_IA_2014=GWAS_IA_2014_multivariate_cov$p_wald,
#                                 p_cov_GA=GWAS_GA_multivariate_cov$p_wald,
#                                 p_cov_stability=GWAS_stability_multivariate_cov_temp$p_wald.y)
# 
# CMplot(GWAS_multivariate_cov,type="p",plot.type="c",chr.labels=paste("Chr",c(1:17),sep=""),r=1,
#        threshold=c(0.05/10937,2.26e-08), threshold.col=c("blue","red"),signal.line=1,signal.col=c("red","green"), 
#        outward=F,cir.chr.h=0.5,chr.den.col="black",file="tiff", signal.cex = 0.6, amplify = T,
#        file.output=T,verbose=TRUE,width=20,height=20)

# FAD samples omitted

GWAS_BC_multivariate_no_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_BC_no_FAD.assoc.txt", header = T)
GWAS_GA_multivariate_no_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_GA_no_FAD.assoc.txt", header = T)
GWAS_IA_multivariate_no_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_IA_2010_no_FAD.assoc.txt", header = T)
GWAS_MN_2015_multivariate_no_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_MN_2015_no_FAD.assoc.txt", header = T)
GWAS_MN_2016_day_1_multivariate_no_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_MN_2016_early_no_FAD.assoc.txt", header = T)
GWAS_MN_2016_day_2_multivariate_no_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_MN_2016_late_no_FAD.assoc.txt", header = T)
GWAS_IA_2013_multivariate_no_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_IA_2013_no_FAD.assoc.txt", header = T)
GWAS_IA_2014_multivariate_no_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_IA_2014_no_FAD.assoc.txt", header = T)
# 
# png("Manhattan_plot_GEMMA_all_envs_multivariate_no_FAD.png", res=120, width = 1500, height = 4000)
# 
# layout(matrix(c(1,2,3,4,5,6,7,8), byrow = F),
#        heights =c(6,6,6,6,6,6,6,7))
# 
# manhattan(GWAS_BC_multivariate_no_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("BC"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08), )
# manhattan(GWAS_GA_multivariate_no_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("GA"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_multivariate_no_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2015_multivariate_no_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2015"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2016_day_1_multivariate_no_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_day_1"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2016_day_2_multivariate_no_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_day_2"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_2013_multivariate_no_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA_2013"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_2014_multivariate_no_FAD, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA_2014"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# dev.off()

# match df lengths
GWAS_IA_multivariate_no_FAD_temp <- merge(GWAS_BC_multivariate_no_FAD,GWAS_IA_multivariate_no_FAD, by = c("chr", "ps"), all.x = T, sort = F)
GWAS_GA_multivariate_no_FAD_temp <- merge(GWAS_BC_multivariate_no_FAD,GWAS_GA_multivariate_no_FAD, by = c("chr", "ps"), all.x = T, sort = F)
GWAS_IA_2013_multivariate_no_FAD_temp <- merge(GWAS_BC_multivariate_no_FAD,GWAS_IA_2013_multivariate_no_FAD, by = c("chr", "ps"), all.x = T, sort = F)
GWAS_IA_2014_multivariate_no_FAD_temp <- merge(GWAS_BC_multivariate_no_FAD,GWAS_IA_2014_multivariate_no_FAD, by = c("chr", "ps"), all.x = T, sort = F)
GWAS_MN_2015_multivariate_no_FAD_temp <- merge(GWAS_BC_multivariate_no_FAD,GWAS_MN_2015_multivariate_no_FAD, by = c("chr", "ps"), all.x = T, sort = F)
GWAS_MN_2016_day_1_multivariate_no_FAD_temp <- merge(GWAS_BC_multivariate_no_FAD,GWAS_MN_2016_day_1_multivariate_no_FAD, by = c("chr", "ps"), all.x = T, sort = F)
GWAS_MN_2016_day_2_multivariate_no_FAD_temp <- merge(GWAS_BC_multivariate_no_FAD,GWAS_MN_2016_day_2_multivariate_no_FAD, by = c("chr", "ps"), all.x = T, sort = F)


GWAS_multivariate_no_FAD <- data.frame(name = rownames(GWAS_BC_multivariate_no_FAD), chr=GWAS_BC_multivariate_no_FAD$chr, snp=GWAS_BC_multivariate_no_FAD$ps, 
                                    p_no_FAD_BC=GWAS_BC_multivariate_no_FAD$p_wald, 
                                    p_no_FAD_MN_2015=GWAS_MN_2015_multivariate_no_FAD_temp$p_wald.y,
                                    p_no_FAD_MN_2016_early=GWAS_MN_2016_day_1_multivariate_no_FAD_temp$p_wald.y,
                                    p_no_FAD_MN_2016_late=GWAS_MN_2016_day_2_multivariate_no_FAD_temp$p_wald.y, 
                                    p_no_FAD_IA_2010=GWAS_IA_multivariate_no_FAD_temp$p_wald.y, 
                                    p_no_FAD_IA_2013=GWAS_IA_2013_multivariate_no_FAD_temp$p_wald.y, 
                                    p_no_FAD_IA_2014=GWAS_IA_2014_multivariate_no_FAD_temp$p_wald.y,
                                    p_no_FAD_GA=GWAS_GA_multivariate_no_FAD_temp$p_wald.y)


CMplot(GWAS_multivariate_no_FAD,type="p",plot.type="c",chr.labels=paste("Chr",c(1:17),sep=""),r=1,
       threshold=c(0.05/10937,2.26e-08), threshold.col=c("blue","red"),signal.line=1,signal.col=c("red","green"), 
       outward=F,cir.chr.h=0.5,chr.den.col="black",file="tiff", signal.cex = 0.6, amplify = T,
       file.output=T,verbose=TRUE,width=20,height=20)

# stderrs

GWAS_BC_multivariate_stderr <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_BC_stderr.assoc.txt", header = T)
GWAS_GA_multivariate_stderr <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_GA_stderr.assoc.txt", header = T)
GWAS_MN_2016_day_1_multivariate_stderr <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_MN_2016_early_stderr.assoc.txt", header = T)
GWAS_MN_2016_day_2_multivariate_stderr <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_MN_2016_late_stderr.assoc.txt", header = T)
GWAS_IA_2013_multivariate_stderr <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_IA_2013_stderr.assoc.txt", header = T)
GWAS_IA_2014_multivariate_stderr <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_IA_2014_stderr.assoc.txt", header = T)
# 
# png("Manhattan_plot_GEMMA_all_envs_multivariate_stderr.png", res=120, width = 1500, height = 3000)
# 
# layout(matrix(c(1,2,3,4,5,6), byrow = F),
#        heights =c(6,6,6,6,6,7))
# 
# manhattan(GWAS_BC_multivariate_stderr, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("BC"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08), )
# manhattan(GWAS_GA_multivariate_stderr, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("GA"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2016_day_1_multivariate_stderr, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_day_1"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_MN_2016_day_2_multivariate_stderr, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("MN_2016_day_2"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_2013_multivariate_stderr, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA_2013"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# manhattan(GWAS_IA_2014_multivariate_stderr, chr = "chr", bp = "ps", snp = "rs", p = "p_wald", main = ("IA_2014"),
#           suggestiveline = -log10(0.05/10937), genomewideline = -log10(2.26e-08))
# dev.off()


GWAS_multivariate_stderr <- data.frame(name = rownames(GWAS_BC_multivariate_stderr), chr=GWAS_BC_multivariate_stderr$chr, snp=GWAS_BC_multivariate_stderr$ps, 
                                p_stderr_BC=GWAS_BC_multivariate_stderr$p_wald, 
                                p_stderr_MN_2016_early=GWAS_MN_2016_day_1_multivariate_stderr$p_wald, 
                                p_stderr_MN_2016_late=GWAS_MN_2016_day_2_multivariate_stderr$p_wald,
                                p_stderr_IA_2013=GWAS_IA_2013_multivariate_stderr$p_wald, 
                                p_stderr_IA_2014=GWAS_IA_2014_multivariate_stderr$p_wald,
                                p_stderr_GA=GWAS_GA_multivariate_stderr$p_wald)
                                #p_stability=GWAS_stability_multivariate_temp$p_wald.y)

# circular plot of multivariate stderr GWAS
CMplot(GWAS_multivariate_stderr,type="p",plot.type="c",chr.labels=paste("Chr",c(1:17),sep=""),r=1,
       threshold=c(0.05/10937,2.26e-08), threshold.col=c("blue","red"),signal.line=1,signal.col=c("red","green"), 
       outward=F,cir.chr.h=0.5,chr.den.col="black",file="tiff", signal.cex = 0.6, amplify = T,
       file.output=T,verbose=TRUE,width=20,height=20)


### create a qq plot to compare all sample GWAS with FAD mutant omission ###

# match df lengths
GWAS_IA_2014_multivariate_no_FAD_temp <- merge(GWAS_IA_2014_multivariate, GWAS_IA_2014_multivariate_no_FAD, by = c("chr", "ps"), all.x = T, sort = F)


GWAS_qqplots <- data.frame(name = rownames(GWAS_IA_2014_multivariate), chr=GWAS_IA_2014_multivariate$chr, snp=GWAS_IA_2014_multivariate$ps, 
                           p_IA=GWAS_IA_2014_multivariate$p_wald, 
                           p_IA_no_FAD=GWAS_IA_2014_multivariate_no_FAD_temp$p_wald.y)

CMplot(GWAS_qqplots,plot.type="q",col=c("dodgerblue1", "darkgoldenrod1"),threshold=2.26e-08,
       ylab.pos=2,signal.pch=c(6,4),signal.cex=1.2,signal.col="red",conf.int=TRUE,box=FALSE,multracks=
               TRUE,cex.axis=2,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,width=10,height=10)

#####################################################
### find lowest p-vals per SNP of all environments###
#####################################################

# oil contents, no FAD covariate
oleic_p_vals <- data.frame(BC = GWAS_BC_oleic$p_wald)
oleic_p_vals$GA <- GWAS_GA_oleic$p_wald
#oleic_p_vals$IA <- GWAS_IA_oleic$p_wald # differing row numbers
oleic_p_vals$IA_2013 <- GWAS_IA_2013_oleic$p_wald
oleic_p_vals$IA_2014 <- GWAS_IA_2014_oleic$p_wald
oleic_p_vals$MN_2015 <- GWAS_MN_2015_oleic$p_wald
oleic_p_vals$MN_2016_early <- GWAS_MN_2016_day_1_oleic$p_wald
oleic_p_vals$MN_2016_late <- GWAS_MN_2016_day_2_oleic$p_wald

min_p_val_oleic <- data.frame(chr = GWAS_BC_oleic$chr, ps = GWAS_BC_oleic$ps, p_wald = apply(oleic_p_vals, 1, min))

palm_p_vals <- data.frame(BC = GWAS_BC_palm$p_wald)
palm_p_vals$GA <- GWAS_GA_palm$p_wald
#palm_p_vals$IA <- GWAS_IA_palm$p_wald
palm_p_vals$IA_2013 <- GWAS_IA_2013_palm$p_wald
palm_p_vals$IA_2014 <- GWAS_IA_2014_palm$p_wald
palm_p_vals$MN_2015 <- GWAS_MN_2015_palm$p_wald
palm_p_vals$MN_2016_early <- GWAS_MN_2016_day_1_palm$p_wald
palm_p_vals$MN_2016_late <- GWAS_MN_2016_day_2_palm$p_wald

min_p_val_palm <- data.frame(chr = GWAS_BC_oleic$chr, ps = GWAS_BC_oleic$ps, p_wald = apply(palm_p_vals, 1, min))

stea_p_vals <- data.frame(BC = GWAS_BC_stea$p_wald)
stea_p_vals$GA <- GWAS_GA_stea$p_wald
#stea_p_vals$IA <- GWAS_IA_stea$p_wald
stea_p_vals$IA_2013 <- GWAS_IA_2013_stea$p_wald
stea_p_vals$IA_2014 <- GWAS_IA_2014_stea$p_wald
stea_p_vals$MN_2015 <- GWAS_MN_2015_stea$p_wald
stea_p_vals$MN_2016_early <- GWAS_MN_2016_day_1_stea$p_wald
stea_p_vals$MN_2016_late <- GWAS_MN_2016_day_2_stea$p_wald

min_p_val_stea <- data.frame(chr = GWAS_BC_oleic$chr, ps = GWAS_BC_oleic$ps, p_wald = apply(stea_p_vals, 1, min))

# oil contents, using FAD covariate

oleic_p_vals <- data.frame(BC = GWAS_BC_oleic$p_wald)
oleic_p_vals$GA <- GWAS_GA_oleic$p_wald
#oleic_p_vals$IA <- GWAS_IA_oleic$p_wald # differing row numbers
oleic_p_vals$IA_2013 <- GWAS_IA_2013_oleic$p_wald
oleic_p_vals$IA_2014 <- GWAS_IA_2014_oleic$p_wald
oleic_p_vals$MN_2015 <- GWAS_MN_2015_oleic$p_wald
oleic_p_vals$MN_2016_early <- GWAS_MN_2016_day_1_oleic$p_wald
oleic_p_vals$MN_2016_late <- GWAS_MN_2016_day_2_oleic$p_wald

min_p_val_oleic <- data.frame(chr = GWAS_BC_oleic$chr, ps = GWAS_BC_oleic$ps, p_wald = apply(oleic_p_vals, 1, min))

palm_FAD_p_vals <- data.frame(BC = GWAS_BC_palm_FAD$p_wald)
palm_FAD_p_vals$GA <- GWAS_GA_palm_FAD$p_wald
#palm_FAD_p_vals$IA <- GWAS_IA_palm_FAD$p_wald
palm_FAD_p_vals$IA_2013 <- GWAS_IA_2013_palm_FAD$p_wald
palm_FAD_p_vals$IA_2014 <- GWAS_IA_2014_palm_FAD$p_wald
palm_FAD_p_vals$MN_2015 <- GWAS_MN_2015_palm_FAD$p_wald
palm_FAD_p_vals$MN_2016_early <- GWAS_MN_2016_day_1_palm_FAD$p_wald
palm_FAD_p_vals$MN_2016_late <- GWAS_MN_2016_day_2_palm_FAD$p_wald

min_p_val_palm_FAD <- data.frame(chr = GWAS_BC_oleic$chr, ps = GWAS_BC_oleic$ps, p_wald = apply(palm_FAD_p_vals, 1, min))

stea_FAD_p_vals <- data.frame(BC = GWAS_BC_stea_FAD$p_wald)
stea_FAD_p_vals$GA <- GWAS_GA_stea_FAD$p_wald
#stea_FAD_p_vals$IA <- GWAS_IA_stea_FAD$p_wald
stea_FAD_p_vals$IA_2013 <- GWAS_IA_2013_stea_FAD$p_wald
stea_FAD_p_vals$IA_2014 <- GWAS_IA_2014_stea_FAD$p_wald
stea_FAD_p_vals$MN_2015 <- GWAS_MN_2015_stea_FAD$p_wald
stea_FAD_p_vals$MN_2016_early <- GWAS_MN_2016_day_1_stea_FAD$p_wald
stea_FAD_p_vals$MN_2016_late <- GWAS_MN_2016_day_2_stea_FAD$p_wald

min_p_val_stea_FAD <- data.frame(chr = GWAS_BC_oleic$chr, ps = GWAS_BC_oleic$ps, p_wald = apply(stea_FAD_p_vals, 1, min))


# oil contents, multivariate
multivariate_p_vals <- data.frame(BC = GWAS_BC_multivariate$p_wald)
multivariate_p_vals$GA <- GWAS_GA_multivariate$p_wald
multivariate_p_vals$IA <- GWAS_IA_multivariate_temp$p_wald.y # differing row numbers fixed before
multivariate_p_vals$IA_2013 <- GWAS_IA_2013_multivariate$p_wald
multivariate_p_vals$IA_2014 <- GWAS_IA_2014_multivariate$p_wald
multivariate_p_vals$MN_2015 <- GWAS_MN_2015_multivariate$p_wald
multivariate_p_vals$MN_2016_early <- GWAS_MN_2016_day_1_multivariate$p_wald
#multivariate_p_vals$MN_2016_late <- GWAS_MN_2016_day_2_multivariate$p_wald # not yet rerun

min_p_val_multivariate <- data.frame(chr = GWAS_BC_multivariate$chr, ps = GWAS_BC_multivariate$ps, p_wald = apply(multivariate_p_vals, 1, min))


# oil stability within environments

oleic_stability_p_vals <- data.frame(BC = GWAS_BC_oleic$p_wald)
oleic_stability_p_vals$GA <- GWAS_GA_oleic$p_wald
oleic_stability_p_vals$IA_2013 <- GWAS_IA_2013_oleic$p_wald
oleic_stability_p_vals$IA_2014 <- GWAS_IA_2014_oleic$p_wald
oleic_stability_p_vals$MN_2016_early <- GWAS_MN_2016_day_1_oleic$p_wald
oleic_stability_p_vals$MN_2016_late <- GWAS_MN_2016_day_2_oleic$p_wald

min_p_val_oleic_stability <- data.frame(chr = GWAS_BC_oleic$chr, ps = GWAS_BC_oleic$ps, p_wald = apply(oleic_stability_p_vals, 1, min))

palm_stability_p_vals <- data.frame(BC = GWAS_BC_palm_stability$p_wald)
palm_stability_p_vals$GA <- GWAS_GA_palm_stability$p_wald
palm_stability_p_vals$IA_2013 <- GWAS_IA_2013_palm_stability$p_wald
palm_stability_p_vals$IA_2014 <- GWAS_IA_2014_palm_stability$p_wald
palm_stability_p_vals$MN_2016_early <- GWAS_MN_2016_day_1_palm_stability$p_wald
palm_stability_p_vals$MN_2016_late <- GWAS_MN_2016_day_2_palm_stability$p_wald

min_p_val_palm_stability <- data.frame(chr = GWAS_BC_oleic$chr, ps = GWAS_BC_oleic$ps, p_wald = apply(palm_stability_p_vals, 1, min))

stea_stability_p_vals <- data.frame(BC = GWAS_BC_stea_stability$p_wald)
stea_stability_p_vals$GA <- GWAS_GA_stea_stability$p_wald
stea_stability_p_vals$IA_2013 <- GWAS_IA_2013_stea_stability$p_wald
stea_stability_p_vals$IA_2014 <- GWAS_IA_2014_stea_stability$p_wald
stea_stability_p_vals$MN_2016_early <- GWAS_MN_2016_day_1_stea_stability$p_wald
stea_stability_p_vals$MN_2016_late <- GWAS_MN_2016_day_2_stea_stability$p_wald

min_p_val_stea_stability <- data.frame(chr = GWAS_BC_oleic$chr, ps = GWAS_BC_oleic$ps, p_wald = apply(stea_stability_p_vals, 1, min))


# min p-vals of multivariate FAD omission
min_p_val_multivariate_no_FAD <- data.frame(chr = GWAS_multivariate_no_FAD$chr, ps = GWAS_multivariate_no_FAD$snp, p_wald = apply(GWAS_multivariate_no_FAD[-c(1:3)], 1, min))

# min p-vals of multivariate stderrs
min_p_val_multivariate_stderr <- data.frame(chr = GWAS_multivariate_stderr$chr, ps = GWAS_multivariate_stderr$snp, p_wald = apply(GWAS_multivariate_stderr[-c(1:3)], 1, min))


# save to avoid loading all separate files again
write.csv(min_p_val_oleic, file = "min_p_val_oleic.csv", row.names = F)
write.csv(min_p_val_palm, file = "min_p_val_palm.csv", row.names = F)
write.csv(min_p_val_stea, file = "min_p_val_stea.csv", row.names = F)
write.csv(min_p_val_multivariate, file = "min_p_val_multivariate.csv", row.names = F)
write.csv(min_p_val_multivariate_no_FAD, file = "min_p_val_multivariate_no_FAD.csv", row.names = F)


write.csv(min_p_val_oleic_FAD, file = "min_p_val_oleic_FAD.csv", row.names = F)
write.csv(min_p_val_palm_FAD, file = "min_p_val_palm_FAD.csv", row.names = F)
write.csv(min_p_val_stea, file = "min_p_val_stea_FAD.csv", row.names = F)

write.csv(min_p_val_oleic_stability, file = "min_p_val_oleic_stability.csv", row.names = F)
write.csv(min_p_val_palm_stability, file = "min_p_val_palm_stability.csv", row.names = F)
write.csv(min_p_val_stea_stability, file = "min_p_val_stea_stability.csv", row.names = F)
write.csv(min_p_val_multivariate_stderr, file = "min_p_val_multivariate_stderr.csv", row.names = F)


######################################
### Generate a candidate gene list ###
######################################

# get significant SNPs
min_p_val_oleic <- read.csv("../GEMMA_vcf2gwas_outputs/min_p_val_oleic.csv")
min_p_val_oleic_subset <- min_p_val_oleic[min_p_val_oleic$p_wald < 0.05/10937,] # subset for significant SNPs

min_p_val_multivariate <- read.csv("min_p_val_multivariate.csv")
min_p_val_multivariate_subset <- min_p_val_multivariate[min_p_val_multivariate$p_wald < 0.05/length(min_p_val_multivariate$chr),] # subset for significant SNPs. Stricter threshold due to p-val inflation

min_p_val_multivariate_no_FAD <- read.csv("min_p_val_multivariate_no_FAD.csv")
min_p_val_multivariate_no_FAD_subset <- min_p_val_multivariate_no_FAD[min_p_val_multivariate_no_FAD$p_wald < 0.05/length(min_p_val_multivariate_no_FAD$chr),] # subset for significant SNPs. Stricter threshold for consistency

min_p_val_multivariate_stderr <- read.csv("min_p_val_multivariate_stderr.csv")
min_p_val_multivariate_stderr_subset <- min_p_val_multivariate_stderr[min_p_val_multivariate_stderr$p_wald < 0.05/length(min_p_val_multivariate_stderr$chr),] # subset for significant SNPs. Stricter threshold for consistency

GWAS_stability_multivariate <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_stability.assoc.txt", header = T)
stability_p_val_subset <- GWAS_stability_multivariate[GWAS_stability_multivariate$p_wald < 0.05/length(GWAS_stability_multivariate$chr),]
#stability_cov_p_val_subset <- GWAS_stability_multivariate_cov[GWAS_stability_multivariate_cov$p_wald < 0.05/length(GWAS_stability_multivariate_cov$chr),] # none above threshold

# run candidate pipeline
multivariate_candidates <- candidate_gene_pipeline(min_p_val_multivariate_subset)
multivariate_no_FAD_candidates <- candidate_gene_pipeline(min_p_val_multivariate_no_FAD_subset)
multivariate_stderr_candidates <- candidate_gene_pipeline(min_p_val_multivariate_stderr_subset)
mulativariate_stability_candidates <- candidate_gene_pipeline(stability_p_val_subset)

# read candidate genes to avoid rerunning pipeline
multivariate_candidates <- read.csv(file = "../multivariate_candidates.csv")
multivariate_no_FAD_candidates <- read.csv(file = "../multivariate_no_FAD_candidates.csv")
multivariate_stderr_candidates <- read.csv(file = "../multivariate_stderr_candidates.csv")
mulativariate_stability_candidates <- read.csv(file = "../multivariate_stability_candidates.csv")


# create Venn diagram of candidate genes between tests
data <- list(multivariate_candidates$DESCRIPTION, multivariate_no_FAD_candidates$DESCRIPTION, multivariate_stderr_candidates$DESCRIPTION, mulativariate_stability_candidates$DESCRIPTION)


# colorblind color palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

venn.diagram(
        x = data,
        category.names = c("mean FAs" , "no FAD" , expression(alpha), expression(beta)),
        filename = 'Figure_10.tiff',
        fill = cbbPalette[1:4],
        cex = 1.5,
        cat.cex = 1.5,
        output=T
)

# write results to tables
write.table(multivariate_candidates, file = "multivariate_candidates.csv", col.names = T, sep = ",")
write.table(multivariate_no_FAD_candidates, file = "multivariate_no_FAD_candidates.csv", col.names = T, sep = ",")
write.table(multivariate_stderr_candidates, file = "multivariate_stderr_candidates.csv", col.names = T, sep = ",")
write.table(mulativariate_stability_candidates, file = "multivariate_stability_candidates.csv", col.names = T, sep = ",")


# accumulate results in one dataframe
multivariate_candidates$test <- "mean FAs"
multivariate_no_FAD_candidates$test <- "mean FAs, FAD omission"
multivariate_stderr_candidates$test <- "alpha stability"
mulativariate_stability_candidates$test <- "beta stability"

df_list <- list(multivariate_candidates, multivariate_no_FAD_candidates, multivariate_stderr_candidates, mulativariate_stability_candidates)

# merge all candidate genes, add column with information in which tests each candidate was found
all_candidates <- df_list %>% reduce(full_join, by='DESCRIPTION')
all_candidates <- merge(merge(merge(multivariate_candidates, multivariate_no_FAD_candidates, by = 1:19, all = T), multivariate_stderr_candidates, by = 1:19, all = T), mulativariate_stability_candidates, by = 1:19, all = T)
all_candidates$tests <- apply(all_candidates[, 20:23] , 1 , paste, collapse = ", ") # create test info column
all_candidates$tests <- gsub(", NA", "", all_candidates$tests)
all_candidates$tests <- gsub("NA, ", "", all_candidates$tests) # get rid of NAs, probably not the best way to do this
all_candidates <- all_candidates[order(all_candidates$seqid, all_candidates$start),] # order by position

all_candidates_output <- all_candidates[c(1,2,4,7,8,10,11,12,24)]
write.table(all_candidates_output, file = "multivariate_candidates_mean_no_FAD_stderrs_beta.csv", col.names = T, sep = ",")
all_candidates_output <- read.csv("multivariate_candidates_mean_no_FAD_stderrs_beta.csv")
print(xtable(all_candidates_output[-c(1,2,6,8)], digits = 0), include.rownames = F)



# compare beta stability with and without FAD omission
GWAS_stability_multivariate <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_stability.assoc.txt", header = T)
GWAS_stability_multivariate_no_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_stability_no_FAD.assoc.txt", header = T)

stability_p_val_subset <- GWAS_stability_multivariate[GWAS_stability_multivariate$p_wald < 0.05/10937,]
stability_no_FAD_p_val_subset <- GWAS_stability_multivariate_no_FAD[GWAS_stability_multivariate_no_FAD$p_wald < 0.05/10937,]

setwd("~/Documents/HDD/Sunflower_stability_GWAS/Post-GWAS_analysis/")
mulativariate_stability_candidates <- candidate_gene_pipeline(stability_p_val_subset)
mulativariate_stability_no_FAD_candidates <- candidate_gene_pipeline(stability_no_FAD_p_val_subset)

