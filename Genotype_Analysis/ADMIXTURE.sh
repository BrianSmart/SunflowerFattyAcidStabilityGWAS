#!/bin/bash

# Run ADMIXTURE using one to ten subpopulations

TARGET_FOLDER=./ADMIXTURE_Results
BED_FILE=/home/briansmart/GATK_SAM_population/IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode_win100_step5_vif10.plink.bed

mkdir $TARGET_FOLDER
cd $TARGET_FOLDER

for K in 1 2 3 4 5 6 7 8 9 10; \
do /home/software/admixture_linux-1.3.0/admixture -j12 --cv $BED_FILE $K | tee log${K}.out; done

# display lowest error K for each number of K
echo "Pick lowest error:"
x=$(grep -h CV log*.out)
echo $x



