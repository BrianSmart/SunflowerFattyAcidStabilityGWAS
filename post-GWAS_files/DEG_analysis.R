setwd("~/Documents/HDD/Sunflower_clean_up/post-GWAS_files/")

# find genes differentially expressed in the seed compared to other tissues
library(DESeq2)

# read expression of different tissues, omitting treatment experiments.
expression_matrix <- as.matrix(read.csv("SRP092742.salmon.count.tsv", sep = "\t", row.names = 1)[c(12,28,29,33,37,46,52:56)]) # SF193_DF_ovary_R1 is seed data

# create sample info file such that ovary is compared to the rest of tissues
sample_info <- data.frame(tissue = colnames(expression_matrix), type = factor(c("treatment", rep("control",10)))) 

# create "DESeqDataSet" object
dds <- DESeqDataSetFromMatrix(countData = round(expression_matrix), colData = DataFrame(sample_info), design = ~ type)

# run DESeq
dds <- DESeq(dds)

# get FDR adjusted p-values
results <- results(dds, alpha = 0.1)
summary(results)

# plot results
plotMA(results) #apparently this normal

library("EnhancedVolcano")


# convert to dataframe
result_df <- as.data.frame(results)
hist(result_df$padj)
hist(result_df$log2FoldChange, breaks = 100)

result_df_005 <- na.omit(result_df[result_df$padj < 0.05,]) # subset for significant DEGs
upregulated_genes_005 <- na.omit(result_df_005[result_df_005$log2FoldChange > 0,])

upregulated_genes_log2 <- na.omit(result_df[result_df$log2FoldChange > 2,])

write.csv(upregulated_genes_005, file = "Seed_DEGs_upregulated_005_pval.csv")
write.csv(upregulated_genes_log2, file = "Seed_DEGs_upregulated_log2.csv")
