# --------------------------- #
#
# Species Distribution Modelling
#
# Jenkins and Stevens (2022) in prep
#
# R script purpose:
# Prepare presence data
# 1. Eunicella verrucosa
# 2. Alyconium digitatum
# 3. Stats for manuscript
#
# Author: Tom Jenkins
# Email: tom.l.jenkins@outlook.com
#
# --------------------------- #

# Load libraries
library(tidyverse)
library(readxl)
library(sf)
library(grid)
library(jpeg)
library(raster)
library(tmap)
library(tmaptools)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(ggpubr)
library(ggspatial)


# --------------- #
#
# Pink sea fan ####
#
# --------------- #

# GBIF query search
# https://www.gbif.org/occurrence/search?offset=0&basis_of_record=HUMAN_OBSERVATION&country=GB&country=FR&country=IE&country=GG&country=JE&country=ES&country=IT&country=PT&country=GI&has_coordinate=true&has_geospatial_issue=false&taxon_key=2263967&occurrence_status=present

# Import pink sea fan occurrences 
pink = read_excel("../data/pinkseafan_allrecords.xlsx")
pink

# Remove Isle of Man record
# https://records.nbnatlas.org/occurrences/8c388555-c4de-46dc-aca4-5c048a8363e2
pink = dplyr::filter(pink, !gbif_occurrenceID %in% 281316412)

# Keep only unique Lon and Lat records
pink = distinct(pink, longitude, latitude, .keep_all = TRUE)
pink

# Check for NA values
map_dbl(pink, ~sum(is.na(.)))

# Convert data.frame to a sf object with LAEA CRS
pink_sf = pink %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = 3035)
pink_sf

# Import Europe basemap from rnaturalworld package
europe = ne_countries(continent = "Europe", scale = "large", returnclass = "sf") %>%
  dplyr::select(name) %>% 
  st_make_valid %>% 
  st_transform(crs = 3035)

# Remove points completely within europe polygon (on land)
on_land = st_within(pink_sf, europe) %>% lengths > 0
pink_sf = pink_sf %>% dplyr::filter(!on_land) 
pink_sf

# Crop to upper Bay of Biscay (extent of Strathclyde rasters)
bb = raster::extent(2277837, 4263921, 2500000, 4650179) %>% st_bbox(crs = 3035)
pink_sf = st_crop(pink_sf, bb)

# Extent of all pink sea fan records
pink_ext = st_bbox(pink_sf)

# Plot all presence records
ggplot()+
  geom_sf(data = europe, colour = "black", fill = "grey", size = 0.1)+
  geom_sf(data = pink_sf, colour = "deeppink")+
  coord_sf(xlim = c(pink_ext$xmin, pink_ext$xmax),
           ylim = c(pink_ext$ymin, pink_ext$ymax),
           expand = TRUE)+
  xlab("Longitude")+
  ylab("Latitude")+
  ggtitle("Study area")

# Export data as csv file with geographic coordinates (WGS84)
pink_coords = pink_sf %>% st_transform(crs = 4326) %>% st_coordinates
pink_sf %>% 
  st_drop_geometry %>% 
  mutate(lon = pink_coords[,"X"], lat = pink_coords[,"Y"]) %>% 
  dplyr::select(species, lon, lat) %>% 
  write_csv("../data/pinkseafan_presence_pts.csv")



# --------------- #
#
# Dead man's fingers ####
#
# --------------- #

# GBIF query search
# https://www.gbif.org/occurrence/map?basis_of_record=HUMAN_OBSERVATION&country=GB&country=IE&country=FR&country=SE&country=NO&country=NL&country=GG&country=JE&country=DE&country=ES&country=PT&country=IS&country=FO&country=DK&country=BE&country=GI&has_coordinate=true&has_geospatial_issue=false&taxon_key=5185175&occurrence_status=present

# Import dead man's fingers occurrences 
dead = read_excel("../data/deadmansfingers_allrecords.xlsx")
dead

# Keep only unique Lon and Lat records
dead = distinct(dead, longitude, latitude, .keep_all = TRUE)
dead

# Check for NA values
map_dbl(dead, ~sum(is.na(.)))

# Convert data.frame to a sf object with LAEA CRS
dead_sf = dead %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = 3035)
dead_sf

# Import Europe basemap from rnaturalworld package
europe = ne_countries(continent = "Europe", scale = "large", returnclass = "sf") %>%
  dplyr::select(name) %>% 
  st_make_valid %>% 
  st_transform(crs = 3035)

# Remove points completely within europe polygon (on land)
on_land = st_within(dead_sf, europe) %>% lengths > 0
dead_sf = dead_sf %>% dplyr::filter(!on_land)
dead_sf

# Crop to upper Bay of Biscay (extent of Strathclyde rasters)
bb = raster::extent(2277837, 4263921, 2500000, 4650179) %>% st_bbox(crs = 3035)
dead_sf = st_crop(dead_sf, bb)

# Extent
dead_ext = st_bbox(dead_sf)

# Plot all presence records
ggplot()+
  geom_sf(data = europe, colour = "black", fill = "grey", size = 0.1)+
  geom_sf(data = dead_sf, colour = "blue")+
  coord_sf(xlim = c(dead_ext$xmin, dead_ext$xmax),
           ylim = c(dead_ext$ymin, dead_ext$ymax),
           expand = TRUE)+
  xlab("Longitude")+
  ylab("Latitude")+
  ggtitle("Study area")

# Export data as csv file with geographic coordinates (WGS84)
dead_coords = dead_sf %>% st_transform(crs = 4326) %>% st_coordinates
dead_sf %>% 
  st_drop_geometry %>% 
  mutate(lon = dead_coords[,"X"], lat = dead_coords[,"Y"]) %>% 
  dplyr::select(species, lon, lat) %>% 
  write_csv("../data/deadmansfingers_presence_pts.csv")



# --------------- #
#
# Figure 1 ####
#
# --------------- #

# Import images
jpeg_psf = readJPEG("../images/E_verrucosa_JRS.jfif") %>% 
  rasterGrob(interpolate = TRUE)
jpeg_psf = ggplot()+
  annotation_custom(jpeg_psf, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_blank()
  )
jpeg_dmf = readJPEG("../images/A_digitatum_JRS.jfif") %>% 
  rasterGrob(interpolate = TRUE)
jpeg_dmf = ggplot()+
  annotation_custom(jpeg_dmf, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_blank()
  )

# Import present-day temperature raster
temperature = raster::stack("../data/raster_predictors/env_rasters.tif") %>% 
  raster::subset("Temp_FromKrige_3km_1951_2000")
raster::plot(temperature)

# Import Europe basemap from rnaturalworld package
europe = ne_countries(continent = "Europe", scale = "large", returnclass = "sf") %>%
  dplyr::select(name) %>% 
  st_make_valid %>% 
  st_transform(crs = 3035)

# Crop to upper Bay of Biscay (extent of Strathclyde rasters)
bb = raster::extent(2377837, 4263921, 2500000, 4650179) %>% st_bbox(crs = 3035)
temperature = raster::crop(temperature, bb)
plot(temperature)
europe = st_crop(europe, bb)
plot(europe)

# Function to plot distribution
plot_distribution = function(background, point_data, geo_extent = NULL, col = col){
  ggplot()+
    geom_sf(data = background, colour = "black", fill = "#d9d9d9", size = 0.1)+
    # layer_spatial(data = background)+
    geom_sf(data = point_data, shape = 21, fill = col, size = 1.5, stroke = 0.25, colour = "black")+
    # coord_sf(expand = TRUE)+
    coord_sf(xlim = c(bb$xmin, bb$xmax), ylim = c(bb$ymin, bb$ymax), expand = FALSE)+
    scale_fill_gradientn(colours = heat.colors(30))+
    annotation_north_arrow(data = europe, location = "bl", height = unit(0.4, "cm"), width = unit(0.4, "cm"),
                           pad_y = unit(0.5, "cm"), style = north_arrow_orienteering(text_size = 4))+
    annotation_scale(data = europe, location = "bl", bar_cols = c("black","white"),
                     height = unit(0.15, "cm"), width_hint = 0.15, text_cex = 0.5)+
    xlab("Longitude")+
    ylab("Latitude")+
    theme(
      axis.text = element_text(colour = "black", size = 6),
      axis.title = element_text(colour = "black", size = 7),
      panel.grid = element_line(colour = "white", size = 0.1),
      panel.background = element_rect(fill = "#deebf7"),
      panel.border = element_rect(fill = NA, colour = "black", size = 0.3),
      plot.title = element_text(size = 10, face = "bold")
    )
}

# Pink sea fan distribution
psf_distrib = plot_distribution(
  background = europe, point_data = pink_sf, geo_extent = pink_ext, col = "deeppink")+
  ggtitle(expression(italic("Eunicella verrucosa")))

# Dead man's fingers distribution
dmf_distrib = plot_distribution(
  background = europe, point_data = dead_sf, col = "royalblue")+
  ggtitle(expression(italic("Alcyonium digitatum")))

# Patchwork
library(patchwork)
plt_list = list(
  psf_distrib,
  dmf_distrib +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank()) 
)
wrap_plots(plt_list, ncol = 2, nrow = 1)
# ggsave(filename = "../figures/Figure1.pdf", width = 8, height = 6)
ggsave(filename = "../figures/Figure1.jpeg", width = 8, height = 6, dpi = 900)


# Figure 1
fig1 = ggarrange(psf_range + annotation_custom(grob = ggplotGrob(jpeg_psf),
                                               xmin = 2277837,
                                               xmax = Inf,
                                               ymin = 2500000,
                                               ymax = Inf)+
                   fig1_theme,
                 dmf_range + annotation_custom(grob = ggplotGrob(jpeg_dmf),
                                               xmin = 2277837,
                                               xmax = Inf,
                                               ymin = 2500000,
                                               ymax = Inf)+
                   fig1_theme
)
fig1
# ggsave(plot = fig1, file = "figures/Figure1.png", width = 8, height = 5, dpi = 600)
# ggsave(plot = fig1, file = "figures/Figure1.pdf", width = 8, height = 5)
