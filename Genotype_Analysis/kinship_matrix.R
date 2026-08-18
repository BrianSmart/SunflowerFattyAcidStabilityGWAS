# create kinship matrix for gemma using VanRaden

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

# create dosage matrix from VCF
# SAM085 & SAM272 missing
VCF <- read.vcfR("IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf")

# extract allele dosage from VCF
dosage <- extract.gt(VCF, as.numeric = F)

# check if genotype and phenotype match
oil_samples <- read.csv("../Phenotype_Analysis/oleic_all_locations.csv")[,1]
oil_samples_numbers <- as.numeric(gsub("PPN","",oil_samples))
oil_samples_numbers %in% as.numeric(gsub("SAM","",dimnames(dosage)[[2]],)) 
# SAM272 is present in phenotype, missing in genotype


# replace VCF notation to numbers. Not elegant but works
dosage <- replace(dosage, dosage == "0|0", 0)
dosage <- replace(dosage, dosage == "0|1", 1)
dosage <- replace(dosage, dosage == "1|0", 1)
dosage <- replace(dosage, dosage == "1|1", 2)


dosage <- as.numeric(dosage)
# create correctly formatted matrix
matrix <- t(matrix(unlist(dosage), nrow = length(dosage)/287, ncol = 287))

# run VanRaden's method to generate kinship matrix
VanRaden_Kinship <- Gmatrix(matrix, method="VanRaden", ploidy=2)

# write kinship to csv
write.table(VanRaden_Kinship, file = "VanRaden_Kinship.csv", row.names = F, col.names = F, sep=",")

# extract inbreeding values from matrix
inbreeding_VR <- data.frame(Nr = 1:length(diag(VanRaden_Kinship)))
inbreeding_VR$inbreeding <- diag(VanRaden_Kinship-1) 
inbreeding_VR$sample <- row.names(VanRaden_Kinship)

# plot inbreeding relative to population
ggplot(inbreeding_VR, aes(Nr, inbreeding)) +
  xlab("") + ylab("Inbreeding coefficient")+
  geom_point(shape=1) +
  theme_bw()

# remove diagonal columns
VanRaden_Kinship_sub <- VanRaden_Kinship
VanRaden_Kinship_sub[lower.tri(VanRaden_Kinship, diag = T)] <- NA


# reformat for pairwise relationships
VanRaden_Kinship_melted <- melt(VanRaden_Kinship_sub, id.var = rownames(VanRaden_Kinship_sub)[1])
VanRaden_Kinship_melted <- na.omit(VanRaden_Kinship_melted)

# heatmap of kinship
ggplot(VanRaden_Kinship_melted, aes(as.factor(Var1), Var2, group=Var2)) + geom_tile(aes(fill = value)) +
  #geom_text(aes(fill = data$value, label = round(data$value, 1))) +
  scale_fill_gradient(low = "yellow", high = "#D6604D", space = "Lab")

# all pairwise kinships
plot(VanRaden_Kinship_melted$value)

# investigate Kinship calculated by GEMMA

GEMMA_Kinship <- read.delim("SAM_Kinship1.cXX.txt", header = F)
#rownames(GEMMA_Kinship) <- 1:288
#colnames(GEMMA_Kinship) <- 1:288

GEMMA_Kinship_sub <- GEMMA_Kinship
GEMMA_Kinship_sub[lower.tri(GEMMA_Kinship_sub, diag = T)] <- NA
GEMMA_Kinship_melted <- melt(GEMMA_Kinship_sub, id.var =rownames(GEMMA_Kinship_sub)[1],na.rm = T)


