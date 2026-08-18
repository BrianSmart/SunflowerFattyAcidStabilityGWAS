# analyse FAD mutation subset VCF

library(pegas)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggrepel)
library(gridExtra)
library(dplyr)
library(AGHmatrix)
library(vcfR)
library(cluster)
library(factoextra)
library(dendextend)
library(ggfortify)
library(ggpubr)
library(poppr)
library(circlize)
library(stringr)

setwd("~/Genotype_Analysis")

VCF <- read.vcfR("IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_FAD_mut.vcf")

maf <- maf(VCF)
hist(maf[,4], breaks=20, xlim = c(0,0.5))
new_maf <- matrix(unlist(str_split(rownames(maf),pattern = "_")),ncol = 2, byrow = T)
maf_df <- data.frame(chr = as.numeric(new_maf[,1]), bp = as.numeric(new_maf[,2]))
maf_df$maf <- maf[,4]
maf_df_14 <- maf_df[maf_df$chr == 14,]
plot(maf_df_14$bp, maf_df_14$maf)


# visualize maf vs SNP position. 
#Maybe use binned average instead
png(filename = "maf_vs_position_fad_mutants.png", width = 3000, height = 600, res = 150)
ggplot(subset(maf_df, chr == 14), aes(x = bp, y = maf)) +
  geom_point(color = "grey")+#, mapping = aes(x=start, y=avg)) +
  geom_smooth(color = "blue", se = F, method = "loess")
  #theme_bw()
  #ggtitle("Chr14")
  #theme(text = element_text(size = 20))
dev.off()


pegas_VCF <- read.vcf("IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_FAD_mut.vcf", to = 2500000)
X <- as.loci(pegas_VCF)
y <- loci2genind(X)
