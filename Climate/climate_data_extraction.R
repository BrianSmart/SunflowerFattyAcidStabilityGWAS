# testing the climate data package used by Innes et al. (2022)

# library(devtools)
# devtools::install_github("jimhester/archive")
# devtools::install_github("mirzacengic/climatedata")

# library(climatedata)
# library(archive)
library(tidyverse)
#library(inborutils)
library(sp)
#library(rgdal)
library(raster)

setwd("/home/markus/Documents/HDD/Sunflower/script_copy/Climate/")

# doesn't work, have to manually download :(
# # setup
# # Get models with all 4 RCP scenarios
# models_all_rcp <- check_models() %>%
#   group_by(model) %>%
#   mutate(n = n()) %>%
#   ungroup() %>%
#   filter(n == 4) %>%
#   distinct(model) %>%
#   pull()
# 
# chelsa <- get_chelsa(type = "bioclim", layer = 1:19, period = c("current"))

# # copied from LILE script for troubleshooting
# env_data <- read.csv("LILE_seed_collection_spreadsheet.csv", header=T) %>% 
#   mutate(source=as.factor(source), population=as.factor(population))
# 
# geo_data <- env_data %>% dplyr::select(source,population,Lat,Long,Elev_m) %>%
#   filter(!source %in% c(2,5,22,32,38)) %>% #remove mistaken/duplicate Appar
#   filter(!is.na(Lat) | !is.na(Long)) %>% #keep only pops that have coordinates (missing coords for source 37, and Appar doesn't have coords)
#   mutate(Lat_s=scale(Lat), Long_s=scale(Long), Elev_m_s=scale(Elev_m)) # scale predictors
# 
# coords <- data.frame(Long=geo_data$Long, Lat=geo_data$Lat,
#                      row.names = geo_data$source) %>% na.omit()

### For my purposes, the monthly cmi (climate moisture index) and tas (average temperature) 
### from chelsa may be best
### to-do: get coordinates and growth months for every location to try regression as GWAS input

# creating a rasterstack with testdata
# chelsa <- raster("CHELSA_clt_01_2012_V.2.1.tif")
# summary(chelsa)
# spplot(chelsa)
# 
# chelsas <- brick(chelsa)
# chelsas[[2]] <- raster("CHELSA_bio10_01.tif")
# chelsas[[3]] <- raster("CHELSA_swb_2015_V.2.1.tif")
# 
# # somewhere in Minnesota
# # Long then Lat, order is important!!!
# coord <- data.frame(Long=-111.20202,Lat=42.5165149)
# 
# points <- SpatialPoints(coord, proj4string = chelsa@crs)
# 
# values <- raster::extract(raster(paste("D:/CHELSA_data_Markus/CHELSA_cmi", "01", "1980", "V.2.1.tif", sep = "_")),points)
# 
# temp <- raster("D:/CHELSA_data_Markus/CHELSA_clt_07_2010_V.2.1.tif")

# extract climate data from monthly CHELSA data, using month, year and 
# coordinates in a dataframe with format data.frame(Long=-111.20202,Lat=42.5165149)
get_climates <- function(month, year, location){
  # this is very slow
  chelsa <- brick(raster(paste("CHELSA_clt", month, year, "V.2.1.tif", sep = "_"), ymx = 83.99986)) # cloud area fraction 
  extent(chelsa) <- extent(-180.0001, 179.9999, -90.00014, 83.99986)
  chelsa[[2]] <- resample(raster(paste("CHELSA_cmi", month, year, "V.2.1.tif", sep = "_")), chelsa) # Climate moisture index
  chelsa[[3]] <- resample(raster(paste("CHELSA_hurs", month, year, "V.2.1.tif", sep = "_")), chelsa) # Near-surface relative humidity
  chelsa[[4]] <- resample(raster(paste("CHELSA_pet", month, year, "V.2.1.tif", sep = "_")), chelsa) # Potential evapotranspiration
  chelsa[[5]] <- resample(raster(paste("CHELSA_pr", month, year, "V.2.1.tif", sep = "_")), chelsa) # Precipitation amount
  chelsa[[6]] <- resample(raster(paste("CHELSA_rsds", month, year, "V.2.1.tif", sep = "_")), chelsa) # Surface downwelling shortwave flux in air
  chelsa[[7]] <- resample(raster(paste("CHELSA_sfcWind", month, year, "V.2.1.tif", sep = "_")), chelsa) # Near-surface wind speed
  chelsa[[8]] <- resample(raster(paste("CHELSA_tas", month, year, "V.2.1.tif", sep = "_")), chelsa) # Mean daily air temperature
  chelsa[[9]] <- resample(raster(paste("CHELSA_tasmax", month, year, "V.2.1.tif", sep = "_")), chelsa) # mean daily maximum air temperature
  chelsa[[10]] <- resample(raster(paste("CHELSA_tasmin", month, year, "V.2.1.tif", sep = "_")), chelsa) # mean daily minimum air temperature
  chelsa[[11]] <- resample(raster(paste("CHELSA_vpd", month, year, "V.2.1.tif", sep = "_")), chelsa) # Vapor pressure deficit
  points <- SpatialPoints(location, proj4string = chelsa@crs)
  return(raster::extract(chelsa,points))
}
  # also very slow
  #location <- SpatialPoints(location, proj4string = chelsa@crs)
  #temp <- raster::extract(brick(raster(paste("CHELSA_clt", month, year, "V.2.1.tif", sep = "_"), ymx = 83.99986)), location)
  #temp <- c(temp, raster::extract(brick(raster(paste("CHELSA_cmi", month, year, "V.2.1.tif", sep = "_"))), location))
  #temp <- c(temp, raster::extract(brick(raster(paste("CHELSA_hurs", month, year, "V.2.1.tif", sep = "_"))), location))
  #temp <- c(temp, raster::extract(brick(raster(paste("CHELSA_pet_penman", month, year, "V.2.1.tif", sep = "_"))), location))
  #temp <- c(temp, raster::extract(brick(raster(paste("CHELSA_pr", month, year, "V.2.1.tif", sep = "_"))), location))
  #temp <- c(temp, raster::extract(brick(raster(paste("CHELSA_rsds", year, month, "V.2.1.tif", sep = "_"))), location))
  #temp <- c(temp, raster::extract(brick(raster(paste("CHELSA_sfcWind", month, year, "V.2.1.tif", sep = "_"))), location))
  #temp <- c(temp, raster::extract(brick(raster(paste("CHELSA_tas", month, year, "V2.1.tif", sep = "_"))), location))
  #temp <- c(temp, raster::extract(brick(raster(paste("CHELSA_tasmax", month, year, "V.2.1.tif", sep = "_"))), location))
  #temp <- c(temp, raster::extract(brick(raster(paste("CHELSA_tasmin", month, year, "V.2.1.tif", sep = "_"))), location))
  #temp <- c(temp, raster::extract(brick(raster(paste("CHELSA_vpd", month, year, "V2.1.tif", sep = "_"))), location))
  #return(temp)

# extract climate data from all time/location combinations

BC_06 <- get_climates("06", "2010", data.frame(Long=-123.250343,Lat=49.256134))
BC_07 <- get_climates("07", "2010", data.frame(Long=-123.250343,Lat=49.256134))
BC_08 <- get_climates("08", "2010", data.frame(Long=-123.250343,Lat=49.256134))

GA_04 <- get_climates("04", "2010", data.frame(Long=-83.530134,Lat=33.875019))
GA_05 <- get_climates("05", "2010", data.frame(Long=-83.530134,Lat=33.875019))
GA_06 <- get_climates("06", "2010", data.frame(Long=-83.530134,Lat=33.875019))

IA_2010_05 <- get_climates("05", "2010", data.frame(Long=-93.663085,Lat=42.010266))
IA_2010_06 <- get_climates("06", "2010", data.frame(Long=-93.663085,Lat=42.010266))
IA_2010_07 <- get_climates("07", "2010", data.frame(Long=-93.663085,Lat=42.010266))

IA_2013_06 <- get_climates("06", "2013", data.frame(Long=-93.663085,Lat=42.010266))
IA_2013_07 <- get_climates("07", "2013", data.frame(Long=-93.663085,Lat=42.010266))
IA_2013_08 <- get_climates("08", "2013", data.frame(Long=-93.663085,Lat=42.010266))

IA_2014_06 <- get_climates("06", "2014", data.frame(Long=-93.663085,Lat=42.010266))
IA_2014_07 <- get_climates("07", "2014", data.frame(Long=-93.663085,Lat=42.010266))
IA_2014_08 <- get_climates("08", "2014", data.frame(Long=-93.663085,Lat=42.010266))

MN_2015_06 <- get_climates("06", "2015", data.frame(Long=-96.63,Lat=46.98))
MN_2015_07 <- get_climates("07", "2015", data.frame(Long=-96.63,Lat=46.98))
MN_2015_08 <- get_climates("08", "2015", data.frame(Long=-96.63,Lat=46.98))

MN_2016_early_05 <- get_climates("05", "2016", data.frame(Long=-96.63,Lat=46.98))
MN_2016_early_06 <- get_climates("06", "2016", data.frame(Long=-96.63,Lat=46.98))
MN_2016_early_07 <- get_climates("07", "2016", data.frame(Long=-96.63,Lat=46.98))

MN_2016_late_06 <- MN_2016_early_06
MN_2016_late_07 <- MN_2016_early_07
MN_2016_late_08 <- get_climates("08", "2016", data.frame(Long=-96.63,Lat=46.98))

climates_temp <- cbind(BC_06, BC_07, BC_08, 
                       GA_04, GA_05, GA_06, 
                       IA_2010_05, IA_2010_06, IA_2010_07, 
                       IA_2013_06, IA_2013_07, IA_2013_08, 
                       IA_2014_06, IA_2014_07, IA_2014_08, 
                       MN_2015_06, MN_2015_07, MN_2015_08, 
                       MN_2016_early_05, MN_2016_early_06, MN_2016_early_07, 
                       MN_2016_late_06, MN_2016_late_07, MN_2016_late_08)

climates <- matrix(unlist(climates_temp), nrow = 8, byrow = T) # reorder climate data so that everything belonging to an environment is in one row

rownames(climates) <- c("BC", "GA", "IA_2010", "IA_2013", "IA_2014", "MN_2015", "MN_2016_early", "MN_2016_late")

# generate column names 
variables <- c("cloud_cover", "climate_moisture", "humidity", "evapotranspiration","precipitation", "air_flux", 
               "wind_speed", "mean_temperature", "max_temperature", "min_temperature", "vapor_pressure")
colnames <- NULL
for (j in 1:3) {
  for (i in variables) {
    colnames <- c(colnames, paste(i, "month", j, sep = "_"))
  }
}
#colnames <- c(colnames, "Lat", "Long")
colnames(climates) <- colnames
climates <- as.data.frame(climates)

# add latitude and longitude
climates$Lat <- c(49.256134, 33.875019, 42.010266, 42.010266, 42.010266, 46.98, 46.98, 46.98)
climates$Long <- c(-123.250343, -83.530134, -93.663085, -93.663085, -93.663085, -96.63, -96.63, -96.63)
#BC_all <- cbind(BC_06, BC_07, BC_08)
#rownames(BC_all) <- c("BC_06","BC_07","BC_08")

# colnames(climates) <- c("clt", "cmi", "hurs", "pet_penman","pr", "rsds", "sfcwind", "tas", "tasmax", "tasmin", "vpd")
# rownames(climates) <- c("BC_06","BC_07","BC_08", "GA_04", "GA_05", "GA_06", "IA_2010_05", "IA_2010_06", "IA_2010_07", 
#                         "IA_2014_06", "IA_2014_07", "IA_2014_08", "MN_2015_06", "MN_2015_07", "MN_2015_08", "MN_2016_early_05", 
#                         "MN_2016_early_06", "MN_2016_early_07", "MN_2016_late_06", "MN_2016_late_07", "MN_2016_late_08")

write.csv(climates, file = "climate_data.csv")
