# Not used for the current paper version: Chernova SNP positions were not lifted over correctly!

setwd("~/Documents/HDD/Sunflower/New_Plots")

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

filtered_list <- list()

for (name in df_names) {
  df <- get(name)
  if (all(c("chr", "ps", "p_wald") %in% names(df))) {
    result <- df %>%
      filter(p_wald < 0.05/nrow(GWAS_BC_multivariate)) %>%
      mutate(source_dataframe = name)
    
    if (nrow(result) > 0) {
      filtered_list[[name]] <- result
    }
  }
}


# Combine all filtered results into one data frame
final_df <- do.call(rbind, filtered_list)

# positions lifted over by Kyle
Chernova_oleic <- data_frame(Chr=c(6,
                                   9,
                                   13,
                                   15), 
                             start=c(68027355,
                                     177501237,
                                     122404993,
                                     40641267
                             ), end=c(68974337,
                                      178093050,
                                      122859111,
                                      40676380)
)
Chernova_oleic$start <- Chernova_oleic$start - 10000000
Chernova_oleic$end <- Chernova_oleic$end + 10000000

Chernova_linoleic <- data_frame(Chr=c(3,
                                      5,
                                      11,
                                      11,
                                      11), 
                                start=c(66345047,
                                        38619256,
                                        52961856,
                                        97251753,
                                        100487580
                                        
                                ), end=c(68407129,
                                         39036273,
                                         53378487,
                                         97677683,
                                         100916849
                                ))

Chernova_linoleic$start <- Chernova_linoleic$start - 10000000
Chernova_linoleic$end <- Chernova_linoleic$end + 10000000

library(GenomicRanges)
gr_snps <- GRanges(seqnames = final_df$chr,
                   ranges = IRanges(start = final_df$ps, end = final_df$ps))

gr_Chernova_oleic <- GRanges(seqnames = Chernova_oleic$Chr,
                             ranges = IRanges(start = Chernova_oleic$start, end = Chernova_oleic$end))

gr_Chernova_linoleic <- GRanges(seqnames = Chernova_linoleic$Chr,
                                ranges = IRanges(start = Chernova_linoleic$start, end = Chernova_linoleic$end))

hits_oleic <- findOverlaps(gr_snps, gr_Chernova_oleic)
snps_in_regions_oleic <- final_df[queryHits(hits_oleic), ]

hits_linoleic <- findOverlaps(gr_snps, gr_Chernova_linoleic)
snps_in_regions_linoleic <- final_df[queryHits(hits_linoleic), ]

write.csv(snps_in_regions_oleic, "Chernova_oleic_dist_less_than_10Mb.csv")
write.csv(snps_in_regions_linoleic, "Chernova_linoleic_dist_less_than_10Mb.csv")
