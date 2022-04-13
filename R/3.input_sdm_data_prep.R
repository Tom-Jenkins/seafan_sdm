# --------------------------- #
#
# Species Distribution Modelling
#
# Jenkins and Stevens (2022) in prep
#
# R script purpose:
#
# 1. Prepare input data for SDMs
# 2. Stats for manuscript
#
# Author: Tom Jenkins
# Email: tom.l.jenkins@outlook.com
#
# --------------------------- #

# Load libraries
library(tidyverse)
library(sf)
library(raster)
library(fasterize)
library(dismo)
library(Rmisc)

# Import bathymetry raster
bathymetry = raster::raster("../data/raster_predictors/bathymetry.tif")
bathymetry
plot(bathymetry)

# Create a bathymetry raster of the continental shelf (<=200 metres)
continental_shelf = bathymetry <= -200
con_shelf = raster::mask(bathymetry, continental_shelf, maskvalue = TRUE)
plot(con_shelf)

# Crop to upper Bay of Biscay (extent of Strathclyde rasters)
ext1 = raster::extent(2277837, 4263921, 2500000, 4650179)
con_shelf = raster::crop(con_shelf, ext1) 
plot(con_shelf)



# --------------- #
#
# Pink sea fan presence data ####
#
# --------------- #

# Import pink sea fan presence points
psf_pres = read_csv("../data/pinkseafan_presence_pts.csv")
psf_pres

# Convert to sf object and transform CRS
psf_pres = st_as_sf(psf_pres, coords = c("lon","lat"), crs = 4326) %>% 
  st_transform(crs = 3035)
psf_pres

# Crop to upper Bay of Biscay (extent of Strathclyde rasters)
psf_pres = st_crop(psf_pres, ext1)

# Plot points on raster
plot(con_shelf, col = topo.colors(100))
points(st_coordinates(psf_pres), pch = 21, bg = "deeppink", col = "black")

# Thin presence records by keeping only one record per raster grid cell
psf_cells = raster::extract(bathymetry, st_coordinates(psf_pres), cellnumbers = TRUE)
psf_cellDups = duplicated(psf_cells[, 1])
psf_pres_thin = psf_pres[!psf_cellDups, ]

# Thin presence records using spThin
spthin_psf = tibble(
  species = psf_pres_thin$species,
  lon = st_coordinates(st_transform(psf_pres_thin, 4326))[, "X"],
  lat = st_coordinates(st_transform(psf_pres_thin, 4326))[, "Y"]
)
set.seed(123)
psf_pres_thin = spThin::thin(
  loc.data = spthin_psf,
  long.col = "lon",
  lat.col = "lat",
  spec.col = "species",
  thin.par = 1,
  reps = 1,
  locs.thinned.list.return = TRUE,
  write.files = FALSE,
  write.log.file = FALSE
  ) %>%
  pluck(1) %>%
  as.data.frame %>%
  st_as_sf(coords = c("Longitude","Latitude"), crs = 4326) %>%
  st_transform(crs = 3035)

# Plot points on raster
plot(con_shelf, col = topo.colors(100))
points(st_coordinates(psf_pres_thin), pch = 21, bg = "deeppink", col = "black")

# Depth ranges of points
raster::extract(bathymetry, psf_pres_thin) %>% summary


# --------------- #
#
# Pink sea fan background points ####
#
# --------------- #

# Restrict background sampling to presence max. depth
psf_max_depth = raster::extract(bathymetry, psf_pres_thin) %>% min(na.rm = TRUE)
psf_max_depth = bathymetry <= psf_max_depth
psf_max_depth = raster::mask(bathymetry, psf_max_depth, maskvalue = TRUE)

# Number of background points
n_bgr_psf = 10000

# Create a 50km buffer around points for sampling background points
# https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.04503
psf_buffer = psf_pres_thin %>%
  # buffer in metres
  st_buffer(., 50000) %>% 
  fasterize(., psf_max_depth)
plot(psf_buffer)

# Remove areas in buffer outside of max. depth
psf_buffer = raster::mask(psf_buffer, psf_max_depth)
plot(psf_buffer)

# Random set of background points using first buffer
set.seed(123)
psf_bgr_pts = randomPoints(
  mask = psf_buffer,
  p = st_coordinates(psf_pres_thin),
  n = n_bgr_psf)

# Depth ranges of presence points
raster::extract(bathymetry, psf_bgr_pts) %>% summary

# Plot points on raster
plot(con_shelf, col = topo.colors(100))
points(psf_bgr_pts, pch = 21, bg = "grey", col = "black")

# Merge presence and background points and export data
psf_pts = tibble(
  pres_backgr = c(rep(1, nrow(psf_pres_thin)), rep(0, nrow(psf_bgr_pts))),
  lon = c(st_coordinates(psf_pres_thin)[,"X"], psf_bgr_pts[,"x"]),
  lat = c(st_coordinates(psf_pres_thin)[,"Y"], psf_bgr_pts[,"y"])
  )
          


# --------------- #
#
# Dead man's fingers presence data ####
#
# --------------- #

# Import dead man's fingers presence points
dmf_pres = read_csv("../data/deadmansfingers_presence_pts.csv")
dmf_pres

# Convert to sf object and transform CRS
dmf_pres = st_as_sf(dmf_pres, coords = c("lon","lat"), crs = 4326) %>% 
  st_transform(crs = 3035)
dmf_pres

# Crop to upper Bay of Biscay (extent of Strathclyde rasters)
dmf_pres = st_crop(dmf_pres, ext1)

# Plot points on raster
plot(con_shelf, col = topo.colors(100))
points(st_coordinates(dmf_pres), pch = 21, bg = "lightblue", col = "black")

# Thin presence records by keeping only one record per raster grid cell
dmf_cells = raster::extract(bathymetry, st_coordinates(dmf_pres), cellnumbers = TRUE)
dmf_cellDups = duplicated(dmf_cells[, 1])
dmf_pres_thin = dmf_pres[!dmf_cellDups, ]

# Thin presence records using spThin
spthin_dmf = tibble(
  species = dmf_pres_thin$species,
  lon = st_coordinates(st_transform(dmf_pres_thin, 4326))[, "X"],
  lat = st_coordinates(st_transform(dmf_pres_thin, 4326))[, "Y"]
)
set.seed(123)
dmf_pres_thin = spThin::thin(
  loc.data = spthin_dmf,
  long.col = "lon",
  lat.col = "lat",
  spec.col = "species",
  thin.par = 1,
  reps = 1,
  locs.thinned.list.return = TRUE,
  write.files = FALSE,
  write.log.file = FALSE
  ) %>% 
  pluck(1) %>% 
  as.data.frame %>% 
  st_as_sf(coords = c("Longitude","Latitude"), crs = 4326) %>% 
  st_transform(crs = 3035)

# Plot points on raster
plot(con_shelf, col = topo.colors(100))
points(st_coordinates(dmf_pres_thin), pch = 21, bg = "lightblue", col = "black")

# Depth ranges of presence points
raster::extract(bathymetry, dmf_pres_thin) %>% summary


# --------------- #
#
# Dead man's fingers background points ####
#
# --------------- #

# Restrict background sampling to presence max. depth
dmf_max_depth = raster::extract(bathymetry, dmf_pres_thin) %>% min(na.rm = TRUE)
dmf_max_depth = bathymetry <= dmf_max_depth
dmf_max_depth = raster::mask(bathymetry, dmf_max_depth, maskvalue = TRUE)

# Number of background points
n_bgr_dmf = 60000

# Create a 50km buffer around points for sampling background points
# https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.04503
dmf_buffer = dmf_pres_thin %>%
  # buffer in metres
  st_buffer(., 50000) %>% 
  fasterize(., dmf_max_depth)
plot(dmf_buffer)

# Remove areas in buffer outside of max. depth
dmf_buffer = raster::mask(dmf_buffer, dmf_max_depth)
plot(dmf_buffer)

# Random set of background points for continental shelf
set.seed(123)
dmf_bgr_pts = randomPoints(
  mask = dmf_buffer,
  p = st_coordinates(dmf_pres_thin),
  n = n_bgr_dmf)

# Plot points on raster
plot(con_shelf, col = topo.colors(100))
points(dmf_bgr_pts, pch = 21, bg = "grey", col = "black")

# Depth ranges of points
raster::extract(bathymetry, dmf_bgr_pts) %>% summary

# Merge presence and background points and export data
dmf_pts = tibble(
  pres_backgr = c(rep(1, nrow(dmf_pres_thin)), rep(0, nrow(dmf_bgr_pts))),
  lon = c(st_coordinates(dmf_pres_thin)[,"X"], dmf_bgr_pts[,"x"]),
  lat = c(st_coordinates(dmf_pres_thin)[,"Y"], dmf_bgr_pts[,"y"])
)


# --------------- #
#
# Extract raster data for points ####
#
# --------------- #

# Combine all points into one tibble
all_pts = tibble(
  lon = c(psf_pts$lon, dmf_pts$lon),
  lat = c(psf_pts$lat, dmf_pts$lat)
)

# Convert to SpatialPoints object and set coordinate reference system
all_pts = SpatialPoints(
  all_pts,
  proj4string = CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
all_pts

# Import all rasters
bathymetry = raster("../data/raster_predictors/bathymetry.tif")
slope = raster("../data/raster_predictors/slope.tif")
strathclyde = raster::stack("../data/raster_predictors/strathclyde.tif")
env_ras = raster::stack("../data/raster_predictors/env_rasters.tif")
# biooracle = raster::stack("../data/raster_predictors/biooracle.tif")

# Combine all raster into one stack
rasters = raster::stack(bathymetry, slope, strathclyde, env_ras)
names(rasters)

# Check CRS matches raster CRS
projection(all_pts) == projection(rasters)

# Extract data for each coordinate
store_data = list()
for (i in 1:nlayers(rasters)){
  store_data[[i]] = raster::extract(rasters[[i]], all_pts)
}

# Combine all data into one tibble
sdm_data = tibble(
  species = c(rep("E_verrucosa", nrow(psf_pts)), rep("A_digitatum", nrow(dmf_pts))),
  pres_backgr = c(psf_pts$pres_backgr, dmf_pts$pres_backgr),
  lon = all_pts@coords[, "lon"],
  lat = all_pts@coords[, "lat"],
)

# Add data to data.frame
names(store_data) = names(rasters)
sdm_data = cbind(sdm_data, as.data.frame(store_data))

# Check for NA values
map_dbl(sdm_data, ~sum(is.na(.)))

# Remove NA records
sdm_data = drop_na(sdm_data)

# Plot density plots of presence versus background points
plot_density = function(df = NULL, species_name = ""){
  df_species = dplyr::select(sdm_data, str_subset(names(sdm_data), "2081_2100", negate = TRUE)) %>% 
    dplyr::filter(., species == {species_name})
  dens_df1 = pivot_longer(df_species, cols = 5:ncol(df_species), names_to = "variable") %>%
    dplyr::filter(pres_backgr == 1)
  dens_df0 = pivot_longer(df_species, cols = 5:ncol(df_species), names_to = "variable") %>%
    dplyr::filter(pres_backgr == 0)
  ggplot()+
    geom_density(data = dens_df1, aes(x = value), color = "red")+
    geom_density(data = dens_df0, aes(x = value), color = "blue")+
    facet_wrap(~variable, scales = "free")+
    ggtitle(species_name)
}
plot_density(sdm_data, species = "E_verrucosa")
plot_density(sdm_data, species = "A_digitatum")

# Export data.frame to csv file
table(sdm_data$species, sdm_data$pres_backgr)
write_csv(sdm_data, file = "../data/sdm_input_data.csv")


# --------------- #
#
# Stats for manuscript ####
#
# --------------- #

# Load libraries
library(tidyverse)

# Import presence data
sdm_data = read_csv("../data/sdm_input_data.csv") %>% 
  dplyr::filter(pres_backgr == 1) %>% 
  dplyr::select(str_subset(names(.), "2081_2100", negate = TRUE))
sdm_data

# Species
Species = "E_verrucosa"
Species = "A_digitatum"

# Variable
# bathymetry | slope | OrbitalVelMean | Rock50cm | TidalVelMean
# Arag_FromKrige_3km_1951_2000 | Calc_FromKrige_3km_1951_2000 | 
# Oxy_FromKrige_3km_1951_2000 | Temp_FromKrige_3km_1951_2000
Variable = "Calc_FromKrige_3km_2081_2100"

# Summary
sdm_data %>% 
  dplyr::filter(species == Species) %>%
  pull({{Variable}}) %>% 
  summary


# --------------- #
#
# Figure 2 ####
#
# --------------- #

# Italic species labels
E_ver = expression(italic("Eunicella verrucosa"))
A_dig = expression(italic("Alcyonium digitatum"))

# Variable facet order
facet_order = c("bathymetry","slope","Rock50cm","OrbitalVelMean","TidalVelMean","Temp_FromKrige_3km_1951_2000","Oxy_FromKrige_3km_1951_2000","Arag_FromKrige_3km_1951_2000","Calc_FromKrige_3km_1951_2000")
facet_labels = c("Bathymetry","Slope","Rock cover","Orbital velocity","Tidal velocity","Temperature","Oxygen concentration","Aragonite saturation state","Calcite saturation state")

# Density plot
sdm_data %>% 
  # transform data to long format for plotting
  pivot_longer(names_to = "Variable", values_to = "Values", cols = 5:15) %>%
  dplyr::filter(!Variable %in% c("Epc_FromKrige_3km_1951_2000","ph_FromKrige_3km_1951_2000")) %>% 
  mutate(facet_var = factor(Variable, levels = facet_order, labels = facet_labels)) %>%
  # plot data
  ggplot(data = ., aes(x = Values, group = species, fill = species))+
  geom_density(alpha = 0.5, adjust = 2)+
  coord_cartesian(expand = TRUE)+
  facet_wrap(~facet_var, scales = "free")+
  scale_fill_manual(values = c("deeppink","royalblue"),
                    limits = c("E_verrucosa","A_digitatum"),
                    labels = c(E_ver,A_dig))+
  ylab("Density")+
  theme_bw()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 12),
    panel.grid = element_line(size = 0.1),
    legend.text = element_text(size = 14),
    legend.position = "top",
    legend.title = element_blank()
  )
ggsave("../figures/Figure2.pdf", width = 11, height = 8)

