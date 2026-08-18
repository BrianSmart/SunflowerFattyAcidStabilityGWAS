setwd("~/Documents/HDD/Sunflower_clean_up/New_Plots")

#beta
GWAS_stability_multivariate <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_stability.assoc.txt", header = T)

# mean FA contents
GWAS_BC_multivariate <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_BC.assoc.txt", header = T)
GWAS_GA_multivariate <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_GA.assoc.txt", header = T)
GWAS_IA_multivariate <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_IA_2010.assoc.txt", header = T)
GWAS_MN_2015_multivariate <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_MN_2015.assoc.txt", header = T)
GWAS_MN_2016_day_1_multivariate <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_MN_2016_early.assoc.txt", header = T)
GWAS_MN_2016_day_2_multivariate <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_MN_2016_late.assoc.txt", header = T)
GWAS_IA_2013_multivariate <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_IA_2013.assoc.txt", header = T)
GWAS_IA_2014_multivariate <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_IA_2014.assoc.txt", header = T)

# FAD samples omitted
GWAS_BC_multivariate_no_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_BC_no_FAD.assoc.txt", header = T)
GWAS_GA_multivariate_no_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_GA_no_FAD.assoc.txt", header = T)
GWAS_IA_multivariate_no_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_IA_2010_no_FAD.assoc.txt", header = T)
GWAS_MN_2015_multivariate_no_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_MN_2015_no_FAD.assoc.txt", header = T)
GWAS_MN_2016_day_1_multivariate_no_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_MN_2016_early_no_FAD.assoc.txt", header = T)
GWAS_MN_2016_day_2_multivariate_no_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_MN_2016_late_no_FAD.assoc.txt", header = T)
GWAS_IA_2013_multivariate_no_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_IA_2013_no_FAD.assoc.txt", header = T)
GWAS_IA_2014_multivariate_no_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_IA_2014_no_FAD.assoc.txt", header = T)

# stderrs
GWAS_BC_multivariate_stderr <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_BC_stderr.assoc.txt", header = T)
GWAS_GA_multivariate_stderr <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_GA_stderr.assoc.txt", header = T)
GWAS_MN_2016_day_1_multivariate_stderr <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_MN_2016_early_stderr.assoc.txt", header = T)
GWAS_MN_2016_day_2_multivariate_stderr <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_MN_2016_late_stderr.assoc.txt", header = T)
GWAS_IA_2013_multivariate_stderr <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_IA_2013_stderr.assoc.txt", header = T)
GWAS_IA_2014_multivariate_stderr <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_IA_2014_stderr.assoc.txt", header = T)



# Get names of all data frames in the global environment
df_names <- ls(envir = .GlobalEnv)
df_names <- df_names[sapply(df_names, function(x) is.data.frame(get(x)))]

library(dplyr)
library(duckplyr)
library(purrr)



# Loop through each data frame
for (name in df_names) {
  df <- get(name)
  
  # Proceed only if 'p_wald' column exists
  if ("p_wald" %in% names(df)) {
    # Filter rows where p_wald < 0.05
    filtered_df <- df[df$p_wald < 0.05/10937, ]
    #filtered_df <- filtered_df[filtered_df]
    
    if (nrow(filtered_df) > 0) {
      filtered_df$source_dataframe <- name  # optional: track origin
      filtered_list[[name]] <- filtered_df
    }
  }
}

# Combine all filtered results into one data frame
final_df <- do.call(rbind, filtered_list)

final_df <- final_df %>% distinct(chr, ps, af, allele1, allele0, p_wald, source_dataframe, .keep_all = F)

# final_df <- final_df %>% distinct(chr, rs, ps, n_miss, allele1, allele0, .keep_all = T)
# final_df <- final_df[c("chr", "ps", "allele1", "allele0")]
# colnames(final_df) <- c("CHROM", "POS", "ALT", "REF")
# final_df$CHROM <- stringr::str_pad(final_df$CHROM, 2, pad = "0")
# 
# readr::write_tsv(final_df, "variants.txt", col_names=FALSE)
# 
final_df$chr_pos <- paste(final_df$chr, final_df$ps, sep = "_")

final_df$source_dataframe <- stringr::str_replace(final_df$source_dataframe, "GWAS_", "") 
final_df$source_dataframe <- stringr::str_replace(final_df$source_dataframe, "_multivariate", "") 
final_df$source_dataframe <- stringr::str_replace(final_df$source_dataframe, "_no_FAD", ", FAD omission") 
final_df$source_dataframe <- stringr::str_replace(final_df$source_dataframe, "_no_FAD", ", FAD omission") 
final_df$source_dataframe <- stringr::str_replace(final_df$source_dataframe, "_day_1", " early") 
final_df$source_dataframe <- stringr::str_replace(final_df$source_dataframe, "_day_2", " late") 
final_df$source_dataframe <- stringr::str_replace(final_df$source_dataframe, "_stderr", " Standard Error model") 
final_df$source_dataframe <- stringr::str_replace(final_df$source_dataframe, "_", " ") 


sort(table(final_df$chr_pos), decreasing = T)

final_df_beta <- final_df[final_df$source_dataframe=="stability",]



presence_table <- final_df %>%
  dplyr::select(chr_pos, source_dataframe) %>%
  distinct() %>%
  mutate(present = TRUE) %>%
  tidyr::pivot_wider(
    names_from = source_dataframe,
    values_from = present,
    values_fill = FALSE
  )


presence_table <- presence_table[order(rowSums(presence_table[-1]), decreasing = T),]


snp_info <- final_df[!duplicated(final_df$chr_pos),]
presence_table_merge <- merge(snp_info[c("chr_pos", "chr", "ps", "allele1", "allele0"#, "af"
                                         )], presence_table, by = "chr_pos")

presence_table_merge <- presence_table_merge[order(rowSums(presence_table_merge[6:27]), decreasing = T),]

presence_table_merge[6:27] <- sapply(presence_table_merge[6:27], function(x) ifelse(x, "X", " "))


# subset to SNPs reaching the less stringent threshold in beta stability w/o FAD
GWAS_stability_multivariate_no_FAD <- read.delim("../GEMMA_vcf2gwas_outputs/GEMMA_stability_no_FAD.assoc.txt", header = T)
GWAS_stability_multivariate_no_FAD$chr_pos <- paste0(GWAS_stability_multivariate_no_FAD$chr, "_", GWAS_stability_multivariate_no_FAD$ps)

stability_no_FAD_p_val_subset <- GWAS_stability_multivariate_no_FAD[GWAS_stability_multivariate_no_FAD$p_wald < 0.05/10937,]

presence_table_stab_overlap <- presence_table_merge[presence_table_merge$chr_pos %in% stability_no_FAD_p_val_subset$chr_pos,]

openxlsx::write.xlsx(presence_table_stab_overlap, "stability_no_FAD_relaxed_p_overlap.xlsx")
