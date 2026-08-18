#!/bin/bash

# month and year to download are input in the command line
month=$1
year=$2

# set the base URL for the data files
base_url="https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/monthly"

# create an array of file types to download
file_types=("clt" "cmi" "hurs" "pet" "pr" "rsds" "sfcWind" "tas" "tasmax" "tasmin" "vpd")

# download each file type
for file_type in "${file_types[@]}"; do
  # construct the URL for the file
  #if [ "$file_type" == "pet" ]; then
  #  # pet has the additional word penman in the file name
  #  file_url="$base_url/$file_type/$year/CHELSA_"$file_type"_penman_"$month"_"$year"_V.2.1.tif"

  #elif [ "$file_type" == "rsds" ]; then
  #  # for rsds, month and year are switched
  #  file_url="$base_url/$file_type/CHELSA_$file_type_"$year"_"$month"_V.2.1.tif"
    
  #else
    file_url="$base_url/$file_type/$year/CHELSA_"$file_type"_"$month"_"$year"_V.2.1.tif"

  #fi
  
  #if [ -e "CHELSA_"$file_type"_penman_"$month"_"$year"_V.2.1.tif" ]; then
  #	continue    # the output file already exists, so skip re-creating it
  #fi
  
  #if [ -e "CHELSA_$file_type_"$year"_"$month"_V.2.1.tif" ]; then
  #  	continue    # the output file already exists, so skip re-creating it
  #fi
  
  if [ -e "CHELSA_"$file_type"_"$month"_"$year"_V.2.1.tif" ]; then
    	continue    # the output file already exists, so skip re-creating it
  fi
  
  # download the file
  wget "$file_url"
done



