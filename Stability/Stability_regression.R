# Analysis of trait stability according to Eberhart & Russell (1966)
### to-do: what to do with genotypes with many missing entries? 
### Especially extreme environments missing may skew the results. Wait for new files
### 277 is only present in IA, it's excluded in every case below

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggfortify)
library(cluster)
library(factoextra)
library(dendextend)
library(readxl)
library(dplyr)
library(pheatmap) 
library(vegan)
library(data.table)
library(gtools)
library(cellWise)
library(metan)

setwd("~/Stability")

# create a function to perform the regression and to extract betai. Intended for use with apply
get_betai <- function(phenotype_list, Ij){
  #if (sum(is.na(phenotype_list)) > 0) { # omit if any phenotype is missing
  if (is.na(phenotype_list[3]) | sum(is.na(phenotype_list)) > 3) { # omit only if GA is missing, OR >3 NAs
    return(NA) # this is necessary to not skew results with missing phenotypes
  }
  else {
  sum <- summary(lm(I(phenotype_list[1:length(phenotype_list)-1]-phenotype_list[length(phenotype_list)]) ~ 0 + Ij)) # regression + summary
  return(sum$coefficients[1]) # extraxt coefficient
  }
}

#########################################
### trait stability over environments ###
#########################################

# oleic phenotype
oleic <- read.csv("../Phenotype_Analysis/oleic_all_locations_updated.csv")
rownames(oleic) <- oleic[,1]
oleic <- oleic[-c(1)] # remove line column
#oleic <- oleic[-c(1,3)] # omit IA for now
#oleic <- oleic[rowSums(is.na(oleic)) != ncol(oleic), ] # omit all NA samples
#oleic <- na.omit(oleic) # omit all NAs as they can heavily skew the regression

oleic[271:287,1] <- NA # blank samples 273 - 288 BC, as there were issues mislabeled fatty acids. Enough to do it here to have it omitted by GEMMA


# calculate regression input values
ui <- apply(oleic, 1, mean, na.rm = T) # will serve as intercepts
mean(as.matrix(oleic), na.rm = T) # should be zero, maybe it's not due to the NAs
Ij <- apply(oleic, 2, mean, na.rm = T)-mean(as.matrix(oleic), na.rm = T) # environmental index

# # test regression to receive betai, used as stability input
# summary <- summary(lm(I(unlist(oleic[5,])-ui[5]) ~ 0 + Ij)) # I sets intercept
# summary$coefficients[1]

plot(y=oleic[5,],x=Ij, ylim = c(0,100))

# add ui to the phenotype dataframe
oleic_and_intercept <- oleic
oleic_and_intercept$intercept <- ui


# apply the function to the phenotypes
oleic_stability <- apply(oleic_and_intercept, 1, get_betai, Ij)

oleic_stability_df <- data.frame(Sample = names(oleic_stability), oleic_stability = oleic_stability)


#linoleic phenotype
linoleic <- read.csv("../Phenotype_Analysis/linoleic_all_locations_updated.csv")

rownames(linoleic) <- linoleic[,1]
linoleic <- linoleic[-c(1)] # remove line column
#linoleic <- linoleic[-c(1,3)] # omit IA for now
#linoleic <- na.omit(linoleic[rowSums(is.na(linoleic)) != ncol(linoleic), ])

# calculate regression input values
ui <- apply(linoleic, 1, mean, na.rm = T) # will serve as intercepts
mean(as.matrix(linoleic), na.rm = T) # sould be zero, maybe it's not due to the NAs
Ij <- apply(linoleic, 2, mean, na.rm = T)-mean(as.matrix(linoleic), na.rm = T) # environmental index

# add ui to the phenotype dataframe
linoleic_and_intercept <- linoleic
linoleic_and_intercept$intercept <- ui

# apply the function to the phenotypes
linoleic_stability <- apply(linoleic_and_intercept, 1, get_betai, Ij)

linoleic_stability_df <- data.frame(Sample = names(linoleic_stability), linoleic_stability = linoleic_stability)


#stearic phenotype
stearic <- read.csv("../Phenotype_Analysis/stearic_all_locations_updated.csv")

rownames(stearic) <- stearic[,1]
stearic <- stearic[-c(1)] # remove line column
#stearic <- stearic[-c(1,3)] # omit IA for now
#stearic <- na.omit(stearic[rowSums(is.na(stearic)) != ncol(stearic), ])

# calculate regression input values
ui <- apply(stearic, 1, mean, na.rm = T) # will serve as intercepts
mean(as.matrix(stearic), na.rm = T) # sould be zero, maybe it's not due to the NAs
Ij <- apply(stearic, 2, mean, na.rm = T)-mean(as.matrix(stearic), na.rm = T) # environmental index

# add ui to the phenotype dataframe
stearic_and_intercept <- stearic
stearic_and_intercept$intercept <- ui

# apply the function to the phenotypes
stearic_stability <- apply(stearic_and_intercept, 1, get_betai, Ij)

stearic_stability_df <- data.frame(Sample = names(stearic_stability), stearic_stability = stearic_stability)


#palmitic phenotype
palmitic <- read.csv("../Phenotype_Analysis/palmitic_all_locations_updated.csv")

rownames(palmitic) <- palmitic[,1]
palmitic <- palmitic[-c(1)] # remove line column
#palmitic <- palmitic[-c(1,3)] # omit IA for now
#palmitic <- na.omit(palmitic[rowSums(is.na(palmitic)) != ncol(palmitic), ])

# calculate regression input values
ui <- apply(palmitic, 1, mean, na.rm = T) # will serve as intercepts
mean(as.matrix(palmitic), na.rm = T) # sould be zero, maybe it's not due to the NAs
Ij <- apply(palmitic, 2, mean, na.rm = T)-mean(as.matrix(palmitic), na.rm = T) # environmental index

# add ui to the phenotype dataframe
palmitic_and_intercept <- palmitic
palmitic_and_intercept$intercept <- ui


# apply the function to the phenotypes
palmitic_stability <- apply(palmitic_and_intercept, 1, get_betai, Ij)

palmitic_stability_df <- data.frame(Sample = names(palmitic_stability), palmitic_stability = palmitic_stability)


# this doesn't handle missing entries well
#all_stabilities <- t(rbind(linoleic_stability, oleic_stability, stearic_stability, palmitic_stability))


# fuse all lists into one dataframe
all_stabilities <- data.frame(linoleic_stability=linoleic_stability,
                              oleic_stability=oleic_stability,
                              palmitic_stability=palmitic_stability,
                              stearic_stability=stearic_stability)

# remove samples with all zeroes, as they result from only one measured phenotype
#test <- all_stabilities[all_stabilities[,1] != 0,]

# check normality
ggqqplot(all_stabilities$oleic)
ggdensity(all_stabilities$oleic)
hist(all_stabilities$oleic_stability, breaks = 20)
hist(all_stabilities$palmitic_stability, breaks = 20)


ggdensity(sqrt(max(all_stabilities$oleic+10)+1)- all_stabilities$oleic+10)
shapiro.test(sqrt(max(all_stabilities$oleic+10)+1)- all_stabilities$oleic+10) # not normal
shapiro.test(transfo(all_stabilities$oleic)$Xt)
shapiro.test(all_stabilities$palmitic_stability)
hist(transfo(all_stabilities$oleic)$Xt, breaks = 20)
test <- transfo(all_stabilities)


# rename samples to match VCF
rownames(all_stabilities) <- gsub("PPN","SAM",rownames(all_stabilities))

# normalize output using Yeo-Johnson power transformation
all_stabilities_normalized <- transfo(all_stabilities)[["Xt"]]



write.table(all_stabilities, file = "oil_stabilities.csv",row.names = T, sep = ",", col.names=NA)
write.table(all_stabilities_normalized, file = "oil_stabilities_normalized.csv",row.names = T, sep = ",", col.names=NA)


# test metan package => I think there's something wrong with the way the function works
# 
# oleic_metan <- oleic
# oleic_metan$Line <- as.factor(rownames(oleic))
# oleic_melt <- melt(oleic_metan)
# oleic_melt$rep <- 1
# oleic_melt[1,4] <- 2
# oleic_melt$variable <- as.numeric(oleic_melt$variable)
# oleic_melt$Line <- as.numeric(oleic_melt$Line)
# 
# metan_oleic <- ge_reg(na.omit(oleic_melt), env = variable, gen = Line, resp = value, rep = rep)

###############################
### trait stderrs over envs ###
###############################

# oleic phenotype
oleic_stderrs <- read.csv("../Phenotype_Analysis/oleic_stderrs_all_locations.csv")
rownames(oleic_stderrs) <- oleic_stderrs[,1]
oleic_stderrs <- oleic_stderrs[-c(1)] # remove line column

oleic_stderrs[271:287,1] <- NA # blank samples 273 - 288 BC, as there were issues mislabeled fatty acids. Enough to do it here to have it omitted by GEMMA


# calculate regression input values
ui <- apply(oleic_stderrs, 1, mean, na.rm = T) # will serve as intercepts
mean(as.matrix(oleic_stderrs), na.rm = T) # should be zero, maybe it's not due to the NAs
Ij <- apply(oleic_stderrs, 2, mean, na.rm = T)-mean(as.matrix(oleic_stderrs), na.rm = T) # environmental index

# # test regression to receive betai, used as stability input
# summary <- summary(lm(I(unlist(oleic[5,])-ui[5]) ~ 0 + Ij)) # I sets intercept
# summary$coefficients[1]

plot(y=oleic_stderrs[20,],x=Ij, ylim = c(0,15))

# add ui to the phenotype dataframe
oleic_stderrs_and_intercept <- oleic_stderrs
oleic_stderrs_and_intercept$intercept <- ui


# apply the function to the phenotypes
oleic_stderrs_stability <- apply(oleic_stderrs_and_intercept, 1, get_betai, Ij)

oleic_stderrs_stability_df <- data.frame(Sample = names(oleic_stderrs_stability), oleic_stderrs_stability = oleic_stderrs_stability)


#linoleic phenotype
linoleic_stderrs <- read.csv("../Phenotype_Analysis/linoleic_stderrs_all_locations.csv")

rownames(linoleic_stderrs) <- linoleic_stderrs[,1]
linoleic_stderrs <- linoleic_stderrs[-c(1)] # remove line column

# calculate regression input values
ui <- apply(linoleic_stderrs, 1, mean, na.rm = T) # will serve as intercepts
mean(as.matrix(linoleic_stderrs), na.rm = T) # sould be zero, maybe it's not due to the NAs
Ij <- apply(linoleic_stderrs, 2, mean, na.rm = T)-mean(as.matrix(linoleic_stderrs), na.rm = T) # environmental index

# add ui to the phenotype dataframe
linoleic_stderrs_and_intercept <- linoleic_stderrs
linoleic_stderrs_and_intercept$intercept <- ui

# apply the function to the phenotypes
linoleic_stderrs_stability <- apply(linoleic_stderrs_and_intercept, 1, get_betai, Ij)

linoleic_stderrs_stability_df <- data.frame(Sample = names(linoleic_stderrs_stability), linoleic_stderrs_stability = linoleic_stderrs_stability)


#stearic phenotype
stearic_stderrs <- read.csv("../Phenotype_Analysis/stearic_stderrs_all_locations.csv")

rownames(stearic_stderrs) <- stearic_stderrs[,1]
stearic_stderrs <- stearic_stderrs[-c(1)] # remove line column

# calculate regression input values
ui <- apply(stearic_stderrs, 1, mean, na.rm = T) # will serve as intercepts
mean(as.matrix(stearic_stderrs), na.rm = T) # sould be zero, maybe it's not due to the NAs
Ij <- apply(stearic_stderrs, 2, mean, na.rm = T)-mean(as.matrix(stearic_stderrs), na.rm = T) # environmental index

# add ui to the phenotype dataframe
stearic_stderrs_and_intercept <- stearic_stderrs
stearic_stderrs_and_intercept$intercept <- ui

# apply the function to the phenotypes
stearic_stderrs_stability <- apply(stearic_stderrs_and_intercept, 1, get_betai, Ij)

stearic_stderrs_stability_df <- data.frame(Sample = names(stearic_stderrs_stability), stearic_stderrs_stability = stearic_stderrs_stability)


#palmitic phenotype
palmitic_stderrs <- read.csv("../Phenotype_Analysis/palmitic_stderrs_all_locations.csv")

rownames(palmitic_stderrs) <- palmitic_stderrs[,1]
palmitic_stderrs <- palmitic_stderrs[-c(1)] # remove line column

# calculate regression input values
ui <- apply(palmitic_stderrs, 1, mean, na.rm = T) # will serve as intercepts
mean(as.matrix(palmitic_stderrs), na.rm = T) # sould be zero, maybe it's not due to the NAs
Ij <- apply(palmitic_stderrs, 2, mean, na.rm = T)-mean(as.matrix(palmitic_stderrs), na.rm = T) # environmental index

# add ui to the phenotype dataframe
palmitic_stderrs_and_intercept <- palmitic_stderrs
palmitic_stderrs_and_intercept$intercept <- ui


# apply the function to the phenotypes
palmitic_stderrs_stability <- apply(palmitic_stderrs_and_intercept, 1, get_betai, Ij)

palmitic_stderrs_stability_df <- data.frame(Sample = names(palmitic_stderrs_stability), palmitic_stderrs_stability = palmitic_stderrs_stability)


# this doesn't handle missing entries well
#all_stabilities <- t(rbind(linoleic_stability, oleic_stability, stearic_stability, palmitic_stability))


# fuse all lists into one dataframe
all_stderrs_stabilities <- data.frame(linoleic_stderrs_stability=linoleic_stderrs_stability,
                              oleic_stderrs_stability=oleic_stderrs_stability,
                              palmitic_stderrs_stability=palmitic_stderrs_stability,
                              stearic_stderrs_stability=stearic_stderrs_stability)

# remove samples with all zeroes, as they result from only one measured phenotype
#test <- all_stabilities[all_stabilities[,1] != 0,]

# check normality
ggqqplot(all_stderrs_stabilities$oleic_stderrs_stability)
ggdensity(all_stderrs_stabilities$oleic_stderrs_stability)
hist(all_stderrs_stabilities$oleic_stderrs_stability, breaks = 20)
hist(all_stderrs_stabilities$palmitic_stderrs_stability, breaks = 20)


#ggdensity(sqrt(max(all_stabilities$oleic+10)+1)- all_stabilities$oleic+10)
shapiro.test(sqrt(max(all_stabilities$oleic+10)+1)- all_stabilities$oleic+10) # not normal
shapiro.test(transfo(all_stabilities$oleic)$Xt)
shapiro.test(all_stderrs_stabilities$oleic_stderrs_stability)
hist(transfo(all_stderrs_stabilities$oleic_stderrs_stability)$Xt, breaks = 20)
shapiro.test(transfo(all_stderrs_stabilities$oleic_stderrs_stability)$Xt)


# rename samples to match VCF
rownames(all_stderrs_stabilities) <- gsub("PPN","SAM",rownames(all_stderrs_stabilities))

# normalize output using Yeo-Johnson power transformation
all_stabilities_normalized <- transfo(all_stderrs_stabilities)[["Xt"]]



write.table(all_stderrs_stabilities, file = "oil_stderrs_stabilities.csv",row.names = T, sep = ",", col.names=NA)
write.table(all_stabilities_normalized, file = "oil_stderrs_stabilities_normalized.csv",row.names = T, sep = ",", col.names=NA)

