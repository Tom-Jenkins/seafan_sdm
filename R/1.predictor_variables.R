# --------------------------- #
#
# Species Distribution Modelling
#
# Jenkins and Stevens (2022) in prep
#
# R script purpose:
# 1. Prepare EMODnet bathymetry data
# 2. Prepare raster predictor variables for SDM
#
# Author: Tom Jenkins
# Email: tom.l.jenkins@outlook.com
#
# --------------------------- #

# Load libraries
library(tidyverse)
library(terra)
library(weathermetrics)


# --------------- #
#
# EMODnet bathymetry prep ####
#
# --------------- #

# Import EMODnet bathymetry data
# https://portal.emodnet-bathymetry.eu/
# file_names = list.files(path = "../0.rasters/EMODnet_bathymetry/raw", pattern = ".asc", full.names = TRUE)
# emod_list = list()
# emod_list = lapply(file_names, function(x) rast(x))

# Merge rasters (very long run time)
# emod_merge = do.call(terra::merge, emod_list)

# Remove land from raster
# temp = emod_merge >= 0
# emod_merge2 = terra::mask(emod_merge, temp, maskvalue = 1)

# Export merged raster
# writeRaster(emod_merge2, filename = "../0.rasters/EMODnet_bathymetry/EMODnet_bathymetry.tif", overwrite = TRUE)

# Import EMODnet bathymetry for Europe
# emod_bathy = rast("../0.rasters/EMODnet_bathymetry/EMODnet_bathymetry.tif")

# LAEA projection https://epsg.io/3035
# laea = "+init=epsg:3035"

# Project to LAEA (equal area grid)
# emod_bathy = project(emod_bathy, laea)
# plot(emod_bathy)

# Export projected raster
# writeRaster(emod_bathy, filename = "../0.rasters/EMODnet_bathymetry/EMODnet_bathymetry_LAEA.tif", overwrite = TRUE)


# --------------- #
#
# Static terrain predictor variables ####
#
# --------------- #

# Import EMODnet bathymetry for Europe
emod_bathy = rast("../0.rasters/EMODnet_bathymetry/EMODnet_bathymetry_LAEA.tif")
plot(emod_bathy)

# Calculate slope from bathymetry grid
slope = terrain(emod_bathy, v = "slope", unit = "degrees")
slope
plot(slope)

# Import Strathclyde rasters
file_list = list.files("../0.rasters/Strathclyde_European_Shelf_Data/",
                      pattern = ".asc", all.files = TRUE, full.names = TRUE)
strathclyde = rast(file_list)
strathclyde

# Import Bio-ORACLE rasters
file_list = list.files("../0.rasters/Bio-Oracle/",
                       pattern = ".tif", all.files = TRUE, full.names = TRUE)
oracle = rast(file_list)

# LAEA projection https://epsg.io/3035
laea = "+init=epsg:3035"

# Project to LAEA
strathclyde_laea = project(strathclyde, laea)
oracle_laea = project(oracle, laea)

# Save extent of Strathclyde rasters
(strath_ext = ext(strathclyde_laea))


# --------------- #
#
# Dynamic predictor variables ####
#
# --------------- #

# Import Morato et al. 2020 raster data
# https://onlinelibrary.wiley.com/doi/10.1111/gcb.14996
# https://doi.pangaea.de/10.1594/PANGAEA.911117?format=html#download
file_list = list.files("../0.rasters/Morato2020_env_variables/",
                       pattern = ".tif", all.files = TRUE, full.names = TRUE)
dynamic_env = rast(file_list)
dynamic_env

# Convert temperature rasters from Kelvin to Celsius
values(dynamic_env$Temp_FromKrige_3km_1951_2000) = values(dynamic_env$Temp_FromKrige_3km_1951_2000, mat = FALSE) %>%  convert_temperature("kelvin", "celsius")
values(dynamic_env$Temp_FromKrige_3km_2081_2100) = values(dynamic_env$Temp_FromKrige_3km_2081_2100, mat = FALSE) %>%  convert_temperature("kelvin", "celsius")
dynamic_env[[9:10]]

# LAEA projection https://epsg.io/3035
laea = "+init=epsg:3035"

# Project to LAEA
dynamic_env_laea = project(dynamic_env, laea)
dynamic_env_laea
plot(dynamic_env_laea[[9:10]])


# --------------- #
#
# Resample raster data ####
#
# --------------- #

# Project raster to the same CRS and resolution as dynamic rasters
bathy_resample = terra::resample(emod_bathy, dynamic_env_laea, method = "bilinear")
bathy_resample
plot(bathy_resample)

# Project raster to the same CRS and resolution as dynamic rasters
slope_resample = terra::resample(slope, dynamic_env_laea, method = "bilinear")
slope_resample
plot(slope_resample)

# Project raster to the same CRS and resolution as dynamic rasters
strath_resample = terra::resample(strathclyde_laea, dynamic_env_laea, method = "bilinear")
strath_resample
plot(strath_resample)

# Project raster to the same CRS and resolution as dynamic rasters
oracle_resample = terra::resample(oracle_laea, dynamic_env_laea, method = "bilinear")
oracle_resample
plot(oracle_resample)

# Crop all rasters to strathclyde extent and plot
# crop(bathy_resample, strath_ext) %>% plot
# crop(slope_resample, strath_ext) %>% plot
# crop(dynamic_env_laea, strath_ext) %>% plot
# crop(strath_resample, strath_ext) %>% plot
crop(oracle_resample, strath_ext) %>% plot


# --------------- #
#
# Export rasters as single geotiff ####
#
# --------------- #

# Export rasters
crop(bathy_resample, strath_ext) %>% 
  writeRaster(., filename = "../data/raster_predictors/bathymetry.tif", overwrite = TRUE)
crop(slope_resample, strath_ext) %>% 
  writeRaster(., filename = "../data/raster_predictors/slope.tif", overwrite = TRUE)
crop(dynamic_env_laea, strath_ext) %>% 
  writeRaster(., filename = "../data/raster_predictors/env_rasters.tif", overwrite = TRUE)
crop(strath_resample, strath_ext) %>% 
  writeRaster(., filename = "../data/raster_predictors/strathclyde.tif", overwrite = TRUE)
crop(oracle_resample, strath_ext) %>% 
  writeRaster(., filename = "../data/raster_predictors/biooracle.tif", overwrite = TRUE)
