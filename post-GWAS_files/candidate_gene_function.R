# the candidate gene functions are defined in this script

library(biomaRt) # load first, as it overwrites some dplyr commands otherwise
#library(qqman)
#library(GWASTools)
#library(SeqArray)
library(ape)
library(zoo)
#library(ggplot2)
#library(ggpubr)
library(dplyr)
#library(qqplotr)
library(stringr)
library(naniar)
library(readxl)
#library(arsenal)
#library(CMplot)
library(tidyr)


# function that extracts the genes within a window of every SNP instead
get_proximate_genes <- function(snp, gff = genes, window = 10000){
  candidate_genes_temp <- gff %>% filter(seqid==snp$chr) %>% 
    #filter(abs(start - snp$ps)==min(abs(start - snp$ps))) %>% 
    filter(start < snp$ps + window & start > snp$ps - window)
  #candidate_genes_temp <- rbind(candidate_genes_temp, gff %>% filter(seqid==snp$chr) %>% filter(abs(start - snp$ps)==min(abs(start - snp$ps))))
  return(candidate_genes_temp)
}

# function for whole candidate gene pipeline
candidate_gene_pipeline <- function(GWAS_output, window = 1e6){
  
  wd <- getwd()
  setwd("~/Documents/HDD/Sunflower_clean_up/post-GWAS_files/")
  
  # get genes around significant SNPs from gff 
  gff <- read.gff("HAN412_Eugene_curated_v1_1.gff3") # gff of the reference genome used, from https://sunflowergenome.org/annotations-data/. 
                                                     # Despite its name, it is the correct file
  gff$seqid <- str_sub(gff$seqid, start = -2)
  gff <- gff[,-(6:8)] # remove unneeded columns
  genes <- subset(gff, type == "gene") # get all genes
  genes <- genes[!duplicated(genes),] # avoid duplicates
  genes$seqid <- as.numeric(genes$seqid)
  
  
  
  # candidate_genes_close <- NULL
  # candidate_genes_close <- get_closest_gene(min_p_val_multivariate_subset[5,], gff = genes) # easy way to initialize the dataframe
  # 
  # i <- 1
  # while (i < length(min_p_val_multivariate_subset$chr)) {
  #   candidate_genes_close <- rbind(candidate_genes_close, get_closest_gene(min_p_val_multivariate_subset[i,], genes))
  #   i <- i+1
  # }
  # candidate_genes_close <- candidate_genes_close[!duplicated(candidate_genes_close),] # filter out duplicate genes
  # 
  
  # in HAN412_Eugene_curated_v1_1.gff3, a gene can be associated with multiple HanXRQv2 genes.
  genes <- na.omit(separate(genes, col = attributes, into = c("rest", "XRQ"), sep = "XRQv2_Gene=")) # create column with XRQ genes, remove genes with no XRQ name
  genes <- separate_rows(genes, XRQ, sep = "-") # separate every XRQ entry into one row
  genes$XRQ <- gsub("HanXRQChr", "HanXRQr2_Chr", genes$XRQ) # match gene names to biomart
  

  # extract INSDC protein IDs using biomaRt to compare with MapMan candidate genes later
  mart<- biomaRt::useEnsemblGenomes(biomart = 'plants_mart',dataset = 'hannuus_eg_gene'#, host = "https://plants.ensembl.org"
                             #, mirror = "asia"
  )
  
  all_gene_names <- biomaRt::getBM(attributes = c("protein_id","ensembl_gene_id","embl"),# ,"uniparc" # extract protein ID and ensembl gene ID. embl matches GeneBank_ID!
                                   #filters="refseq_mrna", # you swap out of this filter for whatever your input is
                                   #values=HanXRQv2_genes, # vector of your NMf
                                   mart=mart)
  temp <- listAttributes(mart = mart)
  
  genes <- merge(all_gene_names, genes, by.x = "ensembl_gene_id", by.y = "XRQ")
  
  genes$protein_id <- sub(pattern = "\\.1",replacement = "", genes$protein_id) %>% tolower()
  
  # read annotation file from https:/mapman.gabipd.org/mapmanstore. 
  # It contains predictions of pathway memberships of genes
  MapMan <- na.omit(read.delim("X4.2_helianthus_annuus.txt", header = T, sep = "\t"))
  MapMan$IDENTIFIER <- gsub("'", "", MapMan$IDENTIFIER) # remove quotation marks from protein ID 
  
  genes_plus_mapman <- merge(genes, MapMan, by.x = "protein_id", by.y = "IDENTIFIER")
  genes_annotated <- genes_plus_mapman[genes_plus_mapman$BINCODE !="'35.1'"& genes_plus_mapman$BINCODE != "'35.2'",] # omit genes without annotation
  
 
  # upregulated DEGs computed using list of expressed genes in seeds from https:/doi.org/10.1038/nature22380 (reference genome paper)
  DEGs <- read.csv(file = "Seed_DEGs_upregulated_log2.csv")
  
  genes_DEG <- merge(genes_annotated, DEGs, by.x = "ensembl_gene_id", by.y = "X")
  
  # to-do: change to look in subset list
  candidate_genes <- NULL
  candidate_genes <- get_proximate_genes(GWAS_output[5,], gff = genes_DEG, 10000) # easy way to initialize the dataframe
  
  i <- 1
  while (i < length(GWAS_output$chr)) {
    candidate_genes <- rbind(candidate_genes, get_proximate_genes(GWAS_output[i,], genes_DEG, window))
    i <- i+1
  }
  candidate_genes <- candidate_genes[!duplicated(candidate_genes),] # filter out duplicate genes
  print(length(candidate_genes[,1]))
  
  duplicate_list <- duplicated(candidate_genes[,c(1,10)])
  candidate_genes_unique <- candidate_genes[!duplicate_list,]
  
  setwd(wd)
  
  return(candidate_genes_unique)
}
