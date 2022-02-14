# --------------------------- #
#
# Species Distribution Modelling
#
# Jenkins and Stevens (2022) in prep
#
# R script purpose:
# Figures for manuscript
#
# Author: Tom Jenkins
# Email: tom.l.jenkins@outlook.com
#
# --------------------------- #

# Load libraries
library(tidyverse)
library(raster)
library(tmap)
library(tmaptools)
library(sf)
library(rgeos)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(ggspatial)
library(jpeg)
library(grid)
library(leaflet)
library(leaflet.providers)
library(leaflet.extras)
library(leaflet.extras2)
library(leafem)
library(htmltools)
library(htmlwidgets)
library(viridis)
library(patchwork)
library(glue)

# Species italics label
psf_lab = expression(italic("Eunicella verrucosa"))
dmf_lab = expression(italic("Alcyonium digitatum"))

# Import UK basemap from rnaturalworld package
# uk = ne_countries(continent = "Europe", scale = "large", returnclass = "sf") %>%
#   dplyr::select(name) %>% 
#   st_transform(crs = 4326)
# uk = crop_shape(uk, bb, polygon = TRUE)
# plot(uk)


# --------------- #
#
# Visualise rasters ####
#
# --------------- #

# Import rasters
ras_predictors = raster::stack("../data/ras_predictors.grd")
names(ras_predictors)

# Merge rasters
ras_static = raster::stack(
  raster::raster("../data/raster_predictors/bathymetry.tif"),
  ras_predictors) %>% 
  raster::subset(str_subset(names(.), "Temp|Oxy|Calc", negate = TRUE))
names(ras_static)

# Import dynamic rasters
ras_dynamic = raster::stack("../data/raster_predictors/env_rasters.tif") %>% 
  raster::subset(str_subset(names(.), "Temp|Oxy|Calc|Arag"))
names(ras_dynamic)

# Plot static raster heatmaps
plt_static = tm_shape(ras_static)+
  tm_raster(title = "", style = "cont")+
  tm_facets(free.scales = TRUE)+
  tm_layout(
    aes.palette = list(seq = "-viridis"),
    legend.position = c("right","bottom"),
    legend.height = 0.3, legend.width = 0.3,
    panel.label.size = 0.9)
plt_static
tmap_save(tm = plt_static, filename = "../figures/FigureS1A.png", width = 5.5, height = 4.5, dpi = 1200)

# Plot static raster heatmaps
plt_dynamic_ras = function(ras_stack){
  tm_shape(ras_stack)+
    tm_raster(title = "", style = "cont")+
    tm_facets(free.scales = FALSE)+
    tm_layout(
      panel.label.size = 0.7,
      legend.outside.size = 0.1)
}
dynamic_vars = c("Temp","Oxy","Arag","Calc")
plt_dynamic_ls = lapply(dynamic_vars, function(x)
  ras_dynamic %>%
    raster::subset(str_subset(names(.), x)) %>% 
    plt_dynamic_ras(.)
    )
plt_dynamic1 = tmap_arrange(
  plt_dynamic_ls[[1]],
  plt_dynamic_ls[[2]],
  ncol = 1)
tmap_save(tm = plt_dynamic1, filename = "../figures/FigureS1Bi.png", width = 7, height = 6, dpi = 1200)
plt_dynamic2 = tmap_arrange(
  plt_dynamic_ls[[3]],
  plt_dynamic_ls[[4]],
  ncol = 1)
tmap_save(tm = plt_dynamic2, filename = "../figures/FigureS1Bii.png", width = 7, height = 6, dpi = 1200)


# --------------- #
#
# Visualise variable contribution ####
#
# --------------- #

# Import data sets
(psf_var = read.csv("../data/pinkseafan_varcontrib.csv"))
(dmf_var = read.csv("../data/deadseafan_varcontrib.csv"))

# Create data.frame of variable contribution
var_contrib = rbind(cbind(psf_var, data.frame(species = "Eunicella")),
                    cbind(dmf_var, data.frame(species = "Alcyonium")))
var_contrib

# Convert to log format
var_contrib_long = pivot_longer(var_contrib, cols = 2:3)
var_contrib_long

# Reorder variables
var_order = c("slope","Rock50cm","OrbitalVelMean","TidalVelMean","Temp_FromKrige_3km_1951_2000","Oxy_FromKrige_3km_1951_2000","Calc_FromKrige_3km_1951_2000")
var_names = c("Slope","Rock cover","Orbital velocity","Tidal velocity","Temperature","Oxygen concentration","Calcite saturation state")
var_contrib_long$variable = factor(var_contrib_long$variable, levels = rev(var_order), labels = rev(var_names))

# Plot variable contribution
plt_varcontrib = function(df, Species = ""){
  df %>% 
    dplyr::filter(species == {Species}) %>% 
    ggplot(data = .)+
      # geom_bar(aes(x = variable, y = value, fill = name),
      #          position = "dodge", stat = "identity",
      #          width = 0.5)+
      geom_col(aes(x = variable, y = value, fill = name),
                position = "dodge", width = 0.7)+
      coord_flip(expand = FALSE)+
      scale_y_continuous(limits = c(0,50))+
      scale_fill_manual(values = c("#d9d9d9","#525252"),
                        labels = c("Percent contribution","Permutation importance"))+
      ylab("(%)")+
      theme_bw()+
      theme(axis.title.y = element_blank(),
            axis.title.x = element_text(size = 8),
            axis.text.y = element_text(size = 9),
            legend.title = element_blank(),
            panel.grid = element_line(size = 0.1),
            plot.title = element_text(size = 10))
}

# Eunicella variable contribution
plt_varA = plt_varcontrib(var_contrib_long, Species = "Eunicella")+
  ggtitle(expression(italic("Eunicella verrucosa")))

# Alcyonium variable contribution
plt_varB = plt_varcontrib(var_contrib_long, Species = "Alcyonium")+
  ggtitle(expression(italic("Alcyonium digitatum")))

# Patchwork
plt_varA +
  plt_varB + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave("../figures/Figure3.pdf", width = 7, height = 4)


#--------------#
#
# Visualise present-day habitat suitability ####
#
#--------------#

# Import raster predictions
psf_raster = raster::stack("../data/pinkseafan_ras.grd")
dmf_raster = raster::stack("../data/deadmansfingers_ras.grd")

# Bounding box
bb = raster::extent(2577837, 4263921, 2500000, 4550179)
bb2 = raster::extent(2777837, 4003921, 2600000, 4232179)

# Plot heatmaps using tmap
plt_habsuit = function(ras, bbox = NULL, labs = "", image_path = NULL){
  ras %>% 
    tm_shape(., bbox = bbox)+
    tm_raster(title = "Habitat suitability", palette = "-inferno")+
    tm_compass(type = "arrow", position = c(0.01,0.06), size = 0.5, text.size = 0.5)+
    tm_scale_bar(position = c(0.01,0.01), breaks = c(0,100,200), lwd = 0.5, text.size = 0.4)+
    tm_layout(
      legend.position = c("RIGHT","BOTTOM"),
      legend.width = 0.22,
      legend.height = 0.22,
      legend.text.size = 0.5,
      title = labs,
      title.position = c("LEFT","TOP"),
      title.size = 0.7
    )+
    tm_logo(image_path, position = c("RIGHT","TOP"), height = 1.8)
}
hs1 = plt_habsuit(psf_raster$psf_pred, bbox = bb2, labs = "Present-day", image_path = "../images/E_verrucosa_JRS.png")
hs2 = plt_habsuit(dmf_raster$dmf_pred, bbox = bb, labs = "Present-day", image_path = "../images/A_digitatum_JRS.png")
hs = tmap_arrange(hs1, hs2)
hs
tmap_save(tm = hs, filename = "../figures/Figure4.png", width = 6.5, height = 4, dpi = 1200)


#--------------#
#
# Interactive Version ####
#
#--------------#

# Import presence points
psf_pres = read_csv("../data/pinkseafan_presence_pts.csv") %>% 
  st_as_sf(., coords = c("lon","lat"), crs = 4326)
dmf_pres = read_csv("../data/deadmansfingers_presence_pts.csv") %>% 
  st_as_sf(., coords = c("lon","lat"), crs = 4326)

# Colours for raster
ras_cols = rev(inferno(5))

# Bins
ras_bin = seq(0, 1, by = 0.2)

# Colour bin
col_bin = colorBin(palette = ras_cols, bins = ras_bin, domain = ras_bin, na.color = "transparent")

# Title
map_title = tags$div(
  HTML("<strong>Jenkins & Stevens 2022: Figure 4</strong>")
) 

# Leaflet map
l1 = leaflet() %>%
  setView(lng = -4.1, lat = 55.0, zoom = 5) %>%
  addProviderTiles(providers$Esri.WorldGrayCanvas) %>%
  # Add a minimap
  addMiniMap(
    position = "topright",
    tiles = providers$OpenStreetMap,
    toggleDisplay = TRUE) %>%
  # Reset map to default setting
  addResetMapButton() %>% 
  # Add a scalebar
  addScaleBar(
    position = "bottomright",
    options = scaleBarOptions(imperial = FALSE)
  ) %>%
  # Add measurement tool
  addMeasure(position = "topleft",
             primaryLengthUnit = "meters",
             secondaryLengthUnit = "kilometers",
             primaryAreaUnit = "sqmeters") %>% 
  # Add pink sea fan raster
  addRasterImage(psf_raster$psf_pred, colors = col_bin, opacity = 0.8,
                 group = "Pink sea fan") %>% 
  # Add pink sea fan presence points
  addCircles(data = psf_pres, color = "black", weight = 1, radius = 500,
             fillColor = "deeppink", fillOpacity = 1,
             group = glue("PSF (N={nrow(psf_pres)})")
  ) %>%
  # Add dead man's fingers raster
  addRasterImage(dmf_raster$dmf_pred, colors = col_bin, opacity = 0.8,
                 group = "Dead man's fingers") %>% 
  # Add dead man's fingers presence points
  addCircles(data = dmf_pres, color = "black", weight = 1, radius = 500,
             fillColor = "royalblue", fillOpacity = 1,
             group = glue("DMF (N={nrow(dmf_pres)})")
  ) %>% 
  # Add layers control
  addLayersControl(
    options = layersControlOptions(collapsed = FALSE),
    baseGroups = c("Pink sea fan","Dead man's fingers"),
    overlayGroups = c(glue("PSF (N={nrow(psf_pres)})"), glue("DMF (N={nrow(dmf_pres)})"))
    ) %>%   
  hideGroup(c("Dead man's fingers", glue("PSF (N={nrow(psf_pres)})"), glue("DMF (N={nrow(dmf_pres)})"))) %>%
  # Add legend
  addLegend(title = "Habitat suitability", pal = col_bin, values = ras_bin, opacity = 0.80) %>% 
  # Add permanent title
  addControl(map_title, position = "bottomleft") %>%
  # Base group title
  htmlwidgets::onRender("
        function() {
            $('.leaflet-control-layers-base').prepend('<label style=\"text-align:left\"><strong>Habitat suitability</strong></label>');
        }
    ") %>% 
  # Overlay group title
  htmlwidgets::onRender("
        function() {
            $('.leaflet-control-layers-overlays').prepend('<label style=\"text-align:left\"><strong>Presence observations</strong></label>');
        }
    ")
l1
saveWidget(l1, file = "../figures/Figure4_interactive.html", title = "Figure4")


#--------------#
#
# Visualise future habitat suitability ####
#
#--------------#

# Plot heatmaps using tmap
plt_habsuit2 = function(ras, bbox = NULL, labs = "", image_path = NULL, compass_bar = TRUE, showlegend = TRUE){
  temp_ras = ras %>% 
    tm_shape(., bbox = bbox)+
    tm_raster(title = "Habitat suitability", palette = "-inferno")+
    tm_layout(
      legend.show = showlegend,
      legend.position = c("RIGHT","BOTTOM"),
      legend.width = 0.22,
      legend.height = 0.22,
      title = labs,
      title.position = c("LEFT","TOP"),
      title.size = 0.6
    )+
    tm_logo(image_path, position = c("RIGHT","TOP"), height = 1.2)
  
  if(compass_bar == TRUE){
    temp_ras = temp_ras+ tm_compass(type = "arrow", position = c(0.01,0.08), size = 0.5, text.size = 0.5)
  }
  if(compass_bar == TRUE){
    temp_ras = temp_ras+ tm_scale_bar(position = c(0.01,0.01), breaks = c(0,100,200), lwd = 0.5, text.size = 0.4)
  }
  return(temp_ras)
}

# Pink sea fan present-day
psf_pre = plt_habsuit2(psf_raster$psf_pred, bbox = bb, image_path = "../images/E_verrucosa_JRS.png", showlegend = FALSE)

# Pink sea fan future
psf_fut = plt_habsuit2(psf_raster$psf_pred_future, bbox = bb, image_path = "../images/E_verrucosa_JRS.png", compass_bar = FALSE)

# Dead man's fingers present-day
dmf_pre = plt_habsuit2(dmf_raster$dmf_pred, bbox = bb, image_path = "../images/A_digitatum_JRS.png", showlegend = FALSE)

# Dead man's fingers future
dmf_fut = plt_habsuit2(dmf_raster$dmf_pred_future, bbox = bb, image_path = "../images/A_digitatum_JRS.png", compass_bar = FALSE)

# Arrange plots
hs_fut = tmap_arrange(ncol = 2,
  psf_pre + tm_layout(title = "Present-day"),
  psf_fut + tm_layout(title = "2081-2100 RCP 8.5"), 
  dmf_pre + tm_layout(title = "Present-day"),
  dmf_fut + tm_layout(title = "2081-2100 RCP 8.5")
  )
# hs_fut
tmap_save(tm = hs_fut, filename = "../figures/Figure5.png", width = 5, height = 6, dpi = 1200)

