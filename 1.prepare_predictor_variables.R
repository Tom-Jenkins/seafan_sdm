# --------------------------- #
#
# Species Distribution Modelling
#
# Jenkins and Stevens (2021) in prep
#
# R script purpose:
# 1. Prepare Bio-ORACLE raster data
# 2. Prepare Strathclyde raster data
# 3. Stats for manuscript
#
# Author: Tom Jenkins
# Email: tom.l.jenkins@outlook.com
# Website: https://tomjenkins.netlify.app
#
# --------------------------- #

# Load libraries
library(tidyverse)
library(sdmpredictors)
library(reshape2)
library(raster)


# ================= #
#
# 1. Prepare Bio-ORACLE raster data ####
#
# ================= #

# --------------- #
#
# Explore Bio-ORACLE datasets ####
#
# --------------- #

# Data available from http://www.bio-oracle.org
# A website containing marine data layers for ecological modelling.
# Files are in .asc format.

# List present-day bio-oracle marine data sets
datasets = list_datasets(terrestrial = FALSE, marine = TRUE)
# write_csv(list_layers(datasets), file = "data/biooracle_layers.csv")
# list_layers_future(datasets)
# list_layers_paleo(datasets)

# Create vectors of temperature layers
temp.bottom = c("BO2_tempmax_bdmean","BO2_tempmean_bdmean","BO2_tempmin_bdmean","BO2_temprange_bdmean")
temp.surface = c("BO2_tempmax_ss","BO2_tempmean_ss","BO2_tempmin_ss","BO2_temprange_ss")
temp_present = c(temp.bottom, temp.surface)

# Create vectors of salinity layers
sal.bottom = c("BO2_salinitymean_bdmean","BO2_salinityrange_bdmean")
sal.surface = c("BO2_salinitymean_ss","BO2_salinityrange_ss")
salinity_present = c(sal.bottom, sal.surface)

# Create a vector of other layers
layers_other = c("MS_bathy_5m","MS_biogeo06_bathy_slope_5m","BO_calcite","BO_ph","BO2_ppmean_bdmean")


# --------------- #
#
# Analyse collinearity of layers ####
#
# --------------- #

# Check correlation between temperature layers
layers_correlation(temp_present) %>% round(digits = 2)
layers_correlation(temp_present) %>% plot_correlation

# Select layers of interest that are not correlated 
temp_present = c("BO2_tempmin_bdmean","BO2_temprange_bdmean","BO2_tempmin_ss","BO2_temprange_ss")
layers_correlation(temp_present) %>% round(digits = 2)
layers_correlation(temp_present) %>% plot_correlation

# Check correlation between salinity layers
layers_correlation(salinity_present) %>% round(digits = 2)
layers_correlation(salinity_present) %>% plot_correlation

# Select layers of interest that are not correlated 
salinity_present = c("BO2_salinitymean_bdmean","BO2_salinityrange_bdmean")
layers_correlation(salinity_present) %>% round(digits = 2)
layers_correlation(salinity_present) %>% plot_correlation

# Check correlation between other layers
layers_correlation(layers_other) %>% round(digits = 2)
layers_correlation(layers_other) %>% plot_correlation

# Change bio-oracle layers to version 2.1
temp_present = gsub("BO2","BO21", temp_present)
salinity_present = gsub("BO2","BO21", salinity_present)
layers_other = gsub("BO2","BO21", layers_other)


# --------------- #
#
# Future temperature layers ####
#
# --------------- #

# Create vectors containing future layers of interest (Bio-ORACLE version 2.1)
temp_future = c("BO21_RCP45_2050_tempmin_bdmean","BO21_RCP85_2050_tempmin_bdmean",
                "BO21_RCP45_2100_tempmin_bdmean","BO21_RCP85_2100_tempmin_bdmean",
                "BO21_RCP45_2050_temprange_bdmean","BO21_RCP85_2050_temprange_bdmean",
                "BO21_RCP45_2100_temprange_bdmean","BO21_RCP85_2100_temprange_bdmean",
                "BO21_RCP45_2050_tempmin_ss","BO21_RCP85_2050_tempmin_ss",
                "BO21_RCP45_2100_tempmin_ss","BO21_RCP85_2100_tempmin_ss",
                "BO21_RCP45_2050_temprange_ss","BO21_RCP85_2050_temprange_ss",
                "BO21_RCP45_2100_temprange_ss","BO21_RCP85_2100_temprange_ss"
)


# --------------- #
#
# Download Bio-ORACLE datasets ####
#
# --------------- #

# Combine present-day and future vectors of variables
layers_temp = c(temp_present, temp_future)
layers_salinity = c(salinity_present)

# Download rasters to sdmpredictors/Meta folder and import into R
options(sdmpredictors_datadir = "C:/R-4.1.0/library/sdmpredictors/Meta/")
rasters_temp = load_layers(layers_temp)
rasters_salinity = load_layers(layers_salinity)
rasters_other = load_layers(layers_other)

# Define a bounding box of the study area
xrange = c(-16, 9)
yrange = c(44.5, 62.75)
bb = extent(xrange, yrange)

# Crop rasters to boundary extent
rasters_temp = crop(rasters_temp, bb)
rasters_salinity = crop(rasters_salinity, bb)
rasters_other = crop(rasters_other, bb)


# --------------- #
#
# Plot rasters ####
#
# --------------- #

# Define colour scheme
cols = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))

# Plot mean bottom minimum temperature
raster::subset(rasters_temp, grep("tempmin_bdmean", names(rasters_temp), value = TRUE)) %>% 
  plot(col = cols(100), zlim = c(0,15), axes = FALSE, box = FALSE)

# Plot mean bottom temperature range
raster::subset(rasters_temp, grep("temprange_bdmean", names(rasters_temp), value = TRUE)) %>% 
  plot(col = cols(100), axes = FALSE, box = FALSE)

# Plot sea surface minimum temperature
raster::subset(rasters_temp, grep("tempmin_ss", names(rasters_temp), value = TRUE)) %>% 
  plot(col = cols(100), axes = FALSE, box = FALSE)

# Plot sea surface temperature range
raster::subset(rasters_temp, grep("temprange_ss", names(rasters_temp), value = TRUE)) %>% 
  plot(col = cols(100), axes = FALSE, box = FALSE)

# Plot mean bottom mean salinity 
raster::subset(rasters_salinity, grep("salinitymean_bdmean", names(rasters_salinity), value = TRUE)) %>% 
  plot(col = cols(100), axes = FALSE, box = FALSE)

# Plot bottom mean salinity range
raster::subset(rasters_salinity, grep("salinityrange_bdmean", names(rasters_salinity), value = TRUE)) %>% 
  plot(col = cols(100), axes = FALSE, box = FALSE)

# Plot other rasters
plot(rasters_other)

# Plot bathymetry using spplot() and limit depth to 200 metres
spplot(rasters_other$MS_bathy_5m, zlim = c(-200,0))


# --------------- #
#
# Combine Bio-ORACLE rasters ####
#
# --------------- #

# Remove temprange rasters
rasters_temp = dropLayer(rasters_temp, grep("range", names(rasters_temp)))

# Combine rasters of interest into one raster stack
ras_biooracle = raster::stack(rasters_temp,
                              # rasters_salinity,
                              rasters_other)
nlayers(ras_biooracle)
names(ras_biooracle)


# ================= #
#
# 2. Prepare Strathclyde raster data ####
#
# ================= #

# Import Strathclyde European shelf rasters
# https://pureportal.strath.ac.uk/en/datasets/data-for-a-synthetic-map-of-the-northwest-european-shelf-sediment
file_list = list.files("data/Strathclyde_European_Shelf_Data/", pattern = ".asc",
                       all.files = TRUE, full.names = TRUE)
ras_strathclyde = raster::stack(file_list)
names(ras_strathclyde)
extent(ras_strathclyde)

# Plot rasters
# plot(ras_strathclyde)

# Change CRS to that of the Bio-ORACLE rasters
crs(ras_strathclyde) = crs(ras_biooracle)

# Resample rasters to have the same resolution as the Bio-ORACLE rasters
ras_strathclyde_rs = resample(ras_strathclyde, ras_biooracle, method = "bilinear")
# plot(ras_strathclyde_rs)

# Convert negative values to zero in orbital velocity and rock cover datasets
ras_strathclyde_rs$OrbitalVelMean[ras_strathclyde_rs$OrbitalVelMean < 0] = 0
ras_strathclyde_rs$Rock50cm[ras_strathclyde_rs$Rock50cm < 0] = 0

# See lines 92-93 on 2.prepare_presence_data.R
# raster::extract(ras_strathclyde_rs$Rock50cm, data.frame(lon = -8.59, lat = 54.83931))
# raster::extract(ras_strathclyde_rs$Rock50cm, data.frame(lon = -8.837907, lat = 54.65946))


# --------------- #
#
# Export raster predictor variables ####
#
# --------------- #

# Combine all rasters into one stack
rasters = raster::stack(ras_biooracle, ras_strathclyde_rs)
rasters

# Export rasters as RData object
save(rasters, file = "data/raster_predictors.RData")


# ================= #
#
# 3. Stats for manuscript ####
#
# ================= #

# --------------- #
#
# Present-day values
#
# --------------- #

# Minimum / Maximum / median
rasters$MS_bathy_5m
rasters$MS_bathy_5m %>% values %>% median(na.rm = TRUE)
rasters$MS_biogeo06_bathy_slope_5m 
rasters$MS_biogeo06_bathy_slope_5m %>% values %>% median(na.rm = TRUE)
rasters$Rock50cm 
rasters$Rock50cm %>% values %>% median(na.rm = TRUE)
rasters$BO21_tempmin_bdmean 
rasters$BO21_tempmin_bdmean %>% values %>% median(na.rm = TRUE)
rasters$OrbitalVelMean 
rasters$OrbitalVelMean %>% values %>% median(na.rm = TRUE)
rasters$TidalVelMean 
rasters$TidalVelMean %>% values %>% median(na.rm = TRUE)
rasters$BO_calcite 
rasters$BO_calcite %>% values %>% median(na.rm = TRUE)
rasters$BO_ph 
rasters$BO_ph %>% values %>% median(na.rm = TRUE)

# --------------- #
#
# Present-day versus projected
#
# --------------- #

# Values of present-day sea bottom temperature
sbt_present = rasters$BO21_tempmin_bdmean %>% values %>% as.vector

# Values of projected sea bottom temperature
sbt_2050_rcp45 = rasters$BO21_RCP45_2050_tempmin_bdmean %>% values %>% as.vector
sbt_2100_rcp45 = rasters$BO21_RCP45_2100_tempmin_bdmean %>% values %>% as.vector

# Create tibble of results
# Remove rows of NA values
sbt_tibble = tibble(sbt_present, sbt_2050_rcp45, sbt_2100_rcp45) %>% 
  drop_na

# Check for NA values
map_dbl(sbt_tibble, ~sum(is.na(.)))

# Minimum / Maximum of datasets
range(sbt_tibble$sbt_present)
range(sbt_tibble$sbt_2050_rcp45)
range(sbt_tibble$sbt_2100_rcp45)

# Calculate differences between present-day and future datasets
sbt_tibble = mutate(sbt_tibble,
                    sbt_diff_2050 = sbt_2050_rcp45 - sbt_present,
                    sbt_diff_2100 = sbt_2100_rcp45 - sbt_present)
sbt_tibble

# Extract the mean difference
mean(sbt_tibble$sbt_diff_2050)
mean(sbt_tibble$sbt_diff_2100)

# Extract the smallest and largest differences
range(sbt_tibble$sbt_diff_2050)
range(sbt_tibble$sbt_diff_2100)


