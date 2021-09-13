# --------------------------- #
#
# Species Distribution Modelling
#
# Jenkins and Stevens (2021) in prep
#
# R script purpose:
# Prepare presence and background points
# 1. Eunicella verrucosa
# 2. Alyconium digitatum
# 3. Stats for manuscript
#
# Author: Tom Jenkins
# Email: tom.l.jenkins@outlook.com
# Website: https://tomjenkins.netlify.app
#
# --------------------------- #

# Load libraries
library(tidyverse)
library(raster)
library(sf)
library(dismo)
library(spThin)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)

# Import rasters
load("data/raster_predictors.RData")
rasters

# Bounding box of the study area
bb = extent(rasters$Rock50cm)

# Import UK basemap from rnaturalworld package
uk = ne_countries(continent = "Europe", scale = "large", returnclass = "sf") %>%
  dplyr::select(name)

# Visualise study area
ggplot()+
  geom_sf(data = uk, colour = "black", fill = "grey", size = 0.5)+
  coord_sf(xlim = c(bb@xmin,bb@xmax), ylim = c(bb@ymin,bb@ymax), expand = FALSE)+
  xlab("Longitude")+
  ylab("Latitude")+
  ggtitle("Study area")


# ================= #
#
# 1. Prepare pink sea fan data ####
#
# ================= #

# --------------- #
#
# Presence data download and filtering ####
#
# --------------- #

# Download pink sea fan records within the study area from the GBIF
# gbif_seafan = gbif(genus = "Eunicella", species = "verrucosa", ext = bb)
# write_csv(gbif_seafan, file = "data/gbif_seafan.csv")
gbif_seafan = read.csv("data/gbif_seafan.csv")
names(gbif_seafan)

# Subset columns
cols_to_keep = c("gbifID","family","species","lon","lat","year","country","institutionCode","identificationVerificationStatus")
gbif_seafan_filt = gbif_seafan %>% dplyr::select(all_of(cols_to_keep))

# Check for NA values
map_dbl(gbif_seafan_filt, ~sum(is.na(.)))

# Remove NA records in lon, lat or institution code
gbif_seafan_filt = gbif_seafan_filt %>% filter(!is.na(lon), !is.na(lat), !is.na(institutionCode))
map_dbl(gbif_seafan_filt, ~sum(is.na(.)))

# Remove records prior to 1990
gbif_seafan_filt = gbif_seafan_filt %>% filter(year >= 1990)

# Keep only unique Lon and Lat records
gbif_seafan_filt = gbif_seafan_filt %>% distinct(lon, lat, .keep_all = TRUE)

# Consider removing "unconfirmed" records (if any)
table(gbif_seafan_filt$identificationVerificationStatus)
# gbif_seafan_filt = gbif_seafan_filt %>%
#   dplyr::filter(!identificationVerificationStatus == "Unconfirmed") %>% 
#   dplyr::filter(!identificationVerificationStatus == "Unconfirmed - plausible") %>% 
#   dplyr::filter(!identificationVerificationStatus == "Unconfirmed - not reviewed")

# Edit -8.510712 longitude so that it falls within the Strathclyde raster limit
# gbif_seafan_filt[1512,]$lon = -8.59

# Print summaries
nrow(gbif_seafan_filt)
range(gbif_seafan_filt$year)
table(gbif_seafan_filt$institutionCode)
table(gbif_seafan_filt$country)

# Plot presence records
ggplot()+
  geom_sf(data = uk, colour = "black", fill = "grey", size = 0.5)+
  geom_point(data = gbif_seafan_filt, aes(x = lon, y = lat), colour = "deeppink")+
  coord_sf(xlim = c(bb@xmin,bb@xmax), ylim = c(bb@ymin,bb@ymax), expand = FALSE)+
  xlab("Longitude")+
  ylab("Latitude")+
  ggtitle("Study area")

# Export gbif_seafan_filt data for use in Figures
write_csv(gbif_seafan_filt, file = "data/gbif_seafan_filt.csv")

# Reduce survey bias in presence records using spatial thinning
# set.seed(123)
# seafan_df_thinned = thin(loc.data = gbif_seafan_filt, lat.col = "lat", long.col = "lon",
#                          spec.col = "species", thin.par = 5, rep = 1,
#                          locs.thinned.list.return = TRUE,
#                          write.files = FALSE, write.log.file = FALSE)

# Retain only one presence record per raster cell to mitigate survey bias
# These will be the final presence points used in the SDMs
seafan_cells = raster::extract(x = rasters$Rock50cm,
                               y = dplyr::select(gbif_seafan_filt, lon, lat),
                               cellnumbers = TRUE)
seafan_cellDups = duplicated(seafan_cells[, 1])
seafan_pres_pts = gbif_seafan_filt[!seafan_cellDups, ] %>%
  as.data.frame %>% 
  dplyr::select(lon, lat)
nrow(seafan_pres_pts)

# Replot pruned presence records on bathymetry raster
plot(rasters$MS_bathy_5m, col = topo.colors(100), main = "MS_bathy_5m")
points(seafan_pres_pts$lon, seafan_pres_pts$lat, pch = 21, bg = "deeppink", col = "black")


# --------------- #
#
# Background points ####
#
# --------------- #

# Define northern and eastern range limit of the pink sea fan
# This ensures random background points are not sampled beyond the
# potential area where pink sea fans are thought to occur or disperse to
seafan_ext = extent(c(xmin = bb@xmin,
                      xmax = max(seafan_pres_pts$lon)+1,
                      ymin = bb@ymin,
                      ymax = max(seafan_pres_pts$lat)
                      )
                    )
seafan_ext

# Extract a random set of n background (pseudo-absence) points from raster
set.seed(123)
seafan_bgr_pts = randomPoints(rasters$Rock50cm,
                              p = seafan_pres_pts,
                              ext = seafan_ext,
                              n = 6000) %>% as.data.frame
colnames(seafan_bgr_pts) = c("lon","lat")

# Remove background points in the North Sea
bgr_remove = subset(seafan_bgr_pts, lon > -1.6 & lat > 52)
seafan_bgr_pts = seafan_bgr_pts[ -c(as.numeric(rownames(bgr_remove))), ]

# Plot background points
ggplot()+
  geom_sf(data = uk, colour = "black", fill = "grey", size = 0.5)+
  geom_point(data = seafan_bgr_pts, aes(x = lon, y = lat), colour = "black")+
  coord_sf(xlim = c(bb@xmin,bb@xmax), ylim = c(bb@ymin,bb@ymax), expand = FALSE)+
  xlab("Longitude")+
  ylab("Latitude")+
  ggtitle("Study area")

# Export seafan_bgr_pts background points for use in Figures
write_csv(seafan_bgr_pts, file = "data/seafan_bgr_pts.csv")

# Check for NA values
map_dbl(seafan_bgr_pts, ~sum(is.na(.)))


# --------------- #
#
# Extract raster data for points ####
#
# --------------- #

# Combine presence and background points
seafan_pts = rbind(seafan_pres_pts, seafan_bgr_pts)

# Add column to indicate presence (1) or background (0) points
seafan_pts$occur.backgr = c(rep(1, nrow(seafan_pres_pts)), rep(0, nrow(seafan_bgr_pts)))

# Convert to SpatialPoints object and set coordinate reference system
seafan_pts = SpatialPoints(seafan_pts, proj4string = CRS("+proj=longlat +datum=WGS84"))
seafan_pts

# Check CRS matches raster CRS
projection(seafan_pts) == projection(rasters)

# Create a data.frame to store environmental data
seafan_data = data.frame(
  Lon = seafan_pts@coords[,"lon"],
  Lat = seafan_pts@coords[,"lat"],
  occur.backgr = seafan_pts@coords[,"occur.backgr"]
)
head(seafan_data)

# Extract data for each coordinate
store_data = list()
for (i in 1:nlayers(rasters)){
  store_data[[i]] = raster::extract(rasters[[i]], seafan_pts)
}

# Add data to data.frame
names(store_data) = names(rasters)
seafan_data = cbind(seafan_data, as.data.frame(store_data))

# Check for NA values
map_dbl(seafan_data, ~sum(is.na(.)))

# Remove NA records
seafan_data = seafan_data %>% drop_na

# Print number of occurrance (1) and background (0) points
table(seafan_data$occur.backgr)

# Round values of environmental variables to three decimal places
# seafan_data[-(1:3)] = apply(seafan_data[-(1:3)], MARGIN = 2, FUN = round, digits = 3)
# seafan_data[1:6, 1:6]

# Replot with final presence and background points on rock cover layer
plot(rasters$Rock50cm, col = topo.colors(100))
points(seafan_data$Lon, seafan_data$Lat, pch = 21, bg = seafan_data$occur.backgr, col = "black")

# Export data.frame to csv file
write_csv(seafan_data, file = "data/sdm_seafan_input_data.csv")


# ================= #
#
# 2. Prepare dead man's fingers data ####
#
# ================= #

# --------------- #
#
# Presence data download and filtering ####
#
# --------------- #

# Download dead man's fingers records within the study area from the GBIF
# gbif_alcyon = gbif(genus = "Alcyonium", species = "digitatum", ext = bb)
# write_csv(gbif_alcyon, file = "data/gbif_alcyon.csv")
gbif_alcyon = read.csv("data/gbif_alcyon.csv")
names(gbif_alcyon)

# Subset columns
cols_to_keep = c("gbifID","family","species","lon","lat","year","country","institutionCode","identificationVerificationStatus")
gbif_alcyon_filt = gbif_alcyon %>% dplyr::select(all_of(cols_to_keep))

# Check for NA values
map_dbl(gbif_alcyon_filt, ~sum(is.na(.)))

# Remove NA records in lon, lat or institution code
gbif_alcyon_filt = gbif_alcyon_filt %>% filter(!is.na(lon), !is.na(lat), !is.na(institutionCode))
map_dbl(gbif_alcyon_filt, ~sum(is.na(.)))

# Remove records prior to 1990
gbif_alcyon_filt = gbif_alcyon_filt %>% filter(year >= 1990)

# Keep only unique Lon and Lat records
gbif_alcyon_filt = gbif_alcyon_filt %>% distinct(lon, lat, .keep_all = TRUE)

# Consider removing "unconfirmed" records (if any)
table(gbif_alcyon_filt$identificationVerificationStatus)
# gbif_alcyon_filt = gbif_alcyon_filt %>%
#   dplyr::filter(!identificationVerificationStatus == "Unconfirmed") %>% 
#   dplyr::filter(!identificationVerificationStatus == "Unconfirmed - plausible") %>% 
#   dplyr::filter(!identificationVerificationStatus == "Unconfirmed - not reviewed")

# Print summaries
nrow(gbif_alcyon_filt)
range(gbif_alcyon_filt$year)
table(gbif_alcyon_filt$institutionCode)
table(gbif_alcyon_filt$country)

# Plot presence records
ggplot()+
  geom_sf(data = uk, colour = "black", fill = "grey", size = 0.5)+
  geom_point(data = gbif_alcyon_filt, aes(x = lon, y = lat), colour = "deeppink")+
  coord_sf(xlim = c(bb@xmin,bb@xmax), ylim = c(bb@ymin,bb@ymax), expand = FALSE)+
  xlab("Longitude")+
  ylab("Latitude")+
  ggtitle("Study area")

# Export gbif_alcyon_filt data for use in Figures
write_csv(gbif_alcyon_filt, file = "data/gbif_alcyon_filt.csv")

# Reduce survey bias in presence records using spatial thinning
# set.seed(123)
# alcyon_df_thinned = thin(loc.data = gbif_alcyon_filt, lat.col = "lat", long.col = "lon",
#                          spec.col = "species", thin.par = 5, rep = 1,
#                          locs.thinned.list.return = TRUE,
#                          write.files = FALSE, write.log.file = FALSE)

# Retain only one presence record per raster cell to mitigate survey bias
# These will be the final presence points used in the SDMs
alcyon_cells = raster::extract(x = rasters$Rock50cm,
                               y = dplyr::select(gbif_alcyon_filt, lon, lat),
                               cellnumbers = TRUE)
alcyon_cellDups = duplicated(alcyon_cells[, 1])
alcyon_pres_pts = gbif_alcyon_filt[!alcyon_cellDups, ] %>%
  as.data.frame %>% 
  dplyr::select(lon, lat)
nrow(alcyon_pres_pts)

# Replot pruned presence records on bathymetry raster
plot(rasters$MS_bathy_5m, col = topo.colors(100), main = "MS_bathy_5m")
points(alcyon_pres_pts$lon, alcyon_pres_pts$lat, pch = 21, bg = "blue", col = "black")


# --------------- #
#
# Background points ####
#
# --------------- #

# Define northern and eastern range limit of the pink sea fan
# This ensures random background points are not sampled beyond the
# potential area where pink sea fans are known to occur or disperse to
# *** 
# This is not required for dead man's fingers because it is found 
# throughout the British Isles and the North Sea
# ***

# Extract a random set of n background (pseudo-absence) points from raster
set.seed(123)
alcyon_bgr_pts = randomPoints(rasters$Rock50cm,
                              p = alcyon_pres_pts,
                              n = 6000) %>% as.data.frame
colnames(alcyon_bgr_pts) = c("lon","lat")

# Plot background points
ggplot()+
  geom_sf(data = uk, colour = "black", fill = "grey", size = 0.5)+
  geom_point(data = alcyon_bgr_pts, aes(x = lon, y = lat), colour = "black")+
  coord_sf(xlim = c(bb@xmin,bb@xmax), ylim = c(bb@ymin,bb@ymax), expand = FALSE)+
  xlab("Longitude")+
  ylab("Latitude")+
  ggtitle("Study area")

# Export alcyon_bgr_pts background points for use in Figures
write_csv(alcyon_bgr_pts, file = "data/alcyon_bgr_pts.csv")

# Check for NA values
map_dbl(alcyon_bgr_pts, ~sum(is.na(.)))


# --------------- #
#
# Extract raster data for points ####
#
# --------------- #

# Combine presence and background points
alcyon_pts = rbind(alcyon_pres_pts, alcyon_bgr_pts)

# Add column to indicate presence (1) or background (0) points
alcyon_pts$occur.backgr = c(rep(1, nrow(alcyon_pres_pts)), rep(0, nrow(alcyon_bgr_pts)))

# Convert to SpatialPoints object and set coordinate reference system
alcyon_pts = SpatialPoints(alcyon_pts, proj4string = CRS("+proj=longlat +datum=WGS84"))
alcyon_pts

# Check CRS matches raster CRS
projection(alcyon_pts) == projection(rasters)

# Create a data.frame to store environmental data
alcyon_data = data.frame(
  Lon = alcyon_pts@coords[,"lon"],
  Lat = alcyon_pts@coords[,"lat"],
  occur.backgr = alcyon_pts@coords[,"occur.backgr"]
)
head(alcyon_data)

# Extract data for each coordinate
store_data = list()
for (i in 1:nlayers(rasters)){
  store_data[[i]] = raster::extract(rasters[[i]], alcyon_pts)
}

# Add data to data.frame
names(store_data) = names(rasters)
alcyon_data = cbind(alcyon_data, as.data.frame(store_data))

# Check for NA values
map_dbl(alcyon_data, ~sum(is.na(.)))

# Remove NA records
alcyon_data = alcyon_data %>% drop_na

# Print number of occurrance (1) and background (0) points
table(alcyon_data$occur.backgr)

# Round values of environmental variables to three decimal places
# alcyon_data[-(1:3)] = apply(alcyon_data[-(1:3)], MARGIN = 2, FUN = round, digits = 3)
# alcyon_data[1:6, 1:6]

# Replot with final presence and background points on rock cover layer
plot(rasters$Rock50cm, col = topo.colors(100))
points(alcyon_data$Lon, alcyon_data$Lat, pch = 21, bg = alcyon_data$occur.backgr, col = "black")

# Export data.frame to csv file
write_csv(alcyon_data, file = "data/sdm_alcyon_input_data.csv")


# ================= #
#
# 3. Stats for manuscript ####
#
# ================= #

# Extract deepest and shallowest depth (presence points only)
seafan_data %>%
  dplyr::filter(occur.backgr == 1) %>% 
  summarise(deepest = min(MS_bathy_5m), shallowest = max(MS_bathy_5m))
alcyon_data %>%
  dplyr::filter(occur.backgr == 1) %>% 
  summarise(deepest = min(MS_bathy_5m), shallowest = max(MS_bathy_5m))

# Calculate % presence points at depth greater or equal to 50 metres
(seafan_data %>% 
  dplyr::filter(occur.backgr == 1) %>%
  dplyr::filter(MS_bathy_5m >= -50) %>% 
  nrow) / table(seafan_data$occur.backgr)[[2]] * 100
(alcyon_data %>% 
  dplyr::filter(occur.backgr == 1) %>%
  dplyr::filter(MS_bathy_5m >= -50) %>% 
  nrow) / table(alcyon_data$occur.backgr)[[2]] * 100

# Extract minimum and maximum slope (presence points only)
seafan_data %>%
  dplyr::filter(occur.backgr == 1) %>% 
  summarise(minimum = min(MS_biogeo06_bathy_slope_5m), maximum = max(MS_biogeo06_bathy_slope_5m))
alcyon_data %>%
  dplyr::filter(occur.backgr == 1) %>% 
  summarise(minimum = min(MS_biogeo06_bathy_slope_5m), maximum = max(MS_biogeo06_bathy_slope_5m))

# Extract minimum and maximum rock cover (presence points only)
seafan_data %>%
  dplyr::filter(occur.backgr == 1) %>% 
  summarise(minimum = min(Rock50cm), maximum = max(Rock50cm))
alcyon_data %>%
  dplyr::filter(occur.backgr == 1) %>% 
  summarise(minimum = min(Rock50cm), maximum = max(Rock50cm))

# Minimum temp
filter(seafan_data, occur.backgr == 1) %>% filter(BO21_tempmin_bdmean < 6.7)

# Extract minimum and maximum sbt (presence points only)
seafan_data %>%
  dplyr::filter(occur.backgr == 1) %>% 
  summarise(minimum = min(BO21_tempmin_bdmean), maximum = max(BO21_tempmin_bdmean))
alcyon_data %>%
  dplyr::filter(occur.backgr == 1) %>% 
  summarise(minimum = min(BO21_tempmin_bdmean), maximum = max(BO21_tempmin_bdmean))

# Extract minimum and maximum orbital velocity (presence points only)
seafan_data %>%
  dplyr::filter(occur.backgr == 1) %>% 
  summarise(minimum = min(OrbitalVelMean), maximum = max(OrbitalVelMean))
alcyon_data %>%
  dplyr::filter(occur.backgr == 1) %>% 
  summarise(minimum = min(OrbitalVelMean), maximum = max(OrbitalVelMean))

# Extract minimum and maximum tidal velocity (presence points only)
seafan_data %>%
  dplyr::filter(occur.backgr == 1) %>% 
  summarise(minimum = min(TidalVelMean), maximum = max(TidalVelMean))
alcyon_data %>%
  dplyr::filter(occur.backgr == 1) %>% 
  summarise(minimum = min(TidalVelMean), maximum = max(TidalVelMean))

# Extract minimum and maximum calcite (presence points only)
seafan_data %>%
  dplyr::filter(occur.backgr == 1) %>% 
  summarise(minimum = min(BO_calcite), maximum = max(BO_calcite))
alcyon_data %>%
  dplyr::filter(occur.backgr == 1) %>% 
  summarise(minimum = min(BO_calcite), maximum = max(BO_calcite))

# Extract minimum and maximum pH (presence points only)
seafan_data %>%
  dplyr::filter(occur.backgr == 1) %>% 
  summarise(minimum = min(BO_ph), maximum = max(BO_ph))
alcyon_data %>%
  dplyr::filter(occur.backgr == 1) %>% 
  summarise(minimum = min(BO_ph), maximum = max(BO_ph))
