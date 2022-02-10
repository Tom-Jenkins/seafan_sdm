# --------------------------- #
#
# Species Distribution Modelling
#
# Jenkins and Stevens (2021) in prep
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
library(RColorBrewer)
library(patchwork)

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
load("../data/ras_predictors.RData")
names(ras_predictors)

# Merge rasters
ras_static = raster::stack(
  raster::raster("../data/raster_predictors/bathymetry.tif"),
  ras_predictors) %>% 
  raster::subset(str_subset(names(.), "Temp|Oxy|Calc", negate = TRUE))
names(rasters)

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
ggsave("../figures/Figure3.jpeg", width = 7, height = 4.5, dpi = 1200)



#--------------#
#
# Visualise present-day habitat suitability ####
#
#--------------#

# Import raster predictions
load("../data/pinkseafan_ras.RData")
load("../data/deadmansfingers_ras.RData")

# Plot heatmaps using tmap
tm_shape(psf_raster)+
  tm_raster(title = "", style = "cont", colorNA = "grey90")+
  tm_facets(free.scales = FALSE)+
  tm_layout(
    # aes.palette = list(seq = "magma")
    )




# Create raster stack
ras_fig3 = stack(ras_pred_filt(seafan_pred$seafan_ras_pred),
                 ras_pred_filt(alcyon_pred$alcyon_ras_pred))

# Facet labels
fig3_facet_lab = c(expression(italic("Eunicella verrucosa")*":"~"present-day"), expression(italic("Alcyonium digitatum")*":"~"present-day"))

# Figure 3
tmap_mode("plot")
fig3 = tm_shape(ras_fig3)+ tm_raster(title = "Habitat suitability", palette = "-RdYlBu")+
  tm_shape(uk, bbox = bb)+ tm_polygons(lwd = 0.2, border.col = "black")+
  tm_compass(type = "arrow", position = c(0.01,0.08), size = 0.7)+
  tm_scale_bar(position = c(0.01,0.01), breaks = c(0,100,200), lwd = 0.5)+
  tm_layout(panel.labels = fig3_facet_lab, panel.label.size = 0.8,
            legend.outside = FALSE, legend.height = 0.15, legend.text.size = 0.7, legend.title.size = 1)
fig3
tmap_save(tm = fig3, filename = "figures/Figure3.png", dpi = 600)
tmap_save(tm = fig3, filename = "figures/Figure3.pdf")


#--------------#
#
# Figure 3: Interactive Version ####
#
#--------------#

# Convert rasters to polygons
ras_fig3_seafan_poly = rasterToPolygons(ras_fig3$seafan_ras_pred)
ras_fig3_alcyon_poly = rasterToPolygons(ras_fig3$alcyon_ras_pred)

# Round habitat suitability to 3 decimal places
round_probs = function(x) { sprintf("%.3f", x) }
ras_fig3_seafan_poly$seafan_ras_pred = round_probs(ras_fig3_seafan_poly$seafan_ras_pred)
ras_fig3_alcyon_poly$alcyon_ras_pred = round_probs(ras_fig3_alcyon_poly$alcyon_ras_pred)

# Convert GBIF presence points to sf object and set CRS
gbif_seafan_filt_sf = st_as_sf(gbif_seafan_filt, coords = c("lon","lat"), crs = 4326)
gbif_alcyon_filt_sf = st_as_sf(gbif_alcyon_filt, coords = c("lon","lat"), crs = 4326)

# Colours for raster
display.brewer.pal(n = 5, name = "RdYlBu")
ras_cols = brewer.pal(5, "RdYlBu") %>% rev

# Bins
ras_bin = seq(0.5, 1.0, by = 0.10)

# Colour bin
col_bin = colorBin(palette = ras_cols, bins = ras_bin, domain = ras_bin, na.color = "transparent")
  
# Figure 3 interactive version
l1 = leaflet() %>%
  # Boundary
  # fitBounds(lng1 = bb@xmin, lat1 = bb@ymin, lng2 = bb@xmax, lat2 = bb@ymax) %>%
  # Add mouse coordinates at top of map
  addMouseCoordinates() %>%
  # Add an inset minimap
  addMiniMap(
    position = "topright",
    tiles = providers$OpenStreetMap,
    toggleDisplay = TRUE) %>%
  # Add logo and link at the top of the map
  leafem::addLogo(
    img = "https://tomjenkins.netlify.app/project/pinkseafan/featured_hu6b938d805fc03dfee7e876a850589a8a_2399056_720x0_resize_q90_lanczos.JPG",
    url = "https://tomjenkins.netlify.app/project/pinkseafan/",
    src = "remote",
    position = "topleft", width = 150, height = 105
  ) %>% 
  # Add a scalebar
  addScaleBar(
    position = "bottomright",
    options = scaleBarOptions(imperial = FALSE)
  ) %>%
  # Reset map to default setting
  addResetMapButton() %>% 
  # Add measurement tool
  addMeasure(position = "topleft",
             primaryLengthUnit = "meters",
             secondaryLengthUnit = "kilometers",
             primaryAreaUnit = "sqmeters") %>% 
  # Add basemap
  addProviderTiles(providers$OpenStreetMap) %>%
  # Add seafan polygons
  addPolygons(data = ras_fig3_seafan_poly, group = "Eunicella verrucosa",
              fillColor = ~col_bin(ras_fig3_seafan_poly$seafan_ras_pred %>% as.numeric),
              fillOpacity = 0.80,
              stroke = 1, weight = 0.2, color = "black",
              popup = ~htmlEscape(ras_fig3_seafan_poly$seafan_ras_pred)) %>%
  # Add seafan presence points
  addCircles(data = gbif_seafan_filt_sf, color = "black", weight = 1, radius = 300,
             fillColor = "deeppink", fillOpacity = 0.8,
             group = "E. verrucosa presence points"
  ) %>% 
  # Add alcyon polygons
  addPolygons(data = ras_fig3_alcyon_poly, group = "Alcyonium digitatum",
              fillColor = ~col_bin(ras_fig3_alcyon_poly$alcyon_ras_pred %>% as.numeric),
              fillOpacity = 0.80,
              stroke = 1, weight = 0.2, color = "black",
              popup = ~htmlEscape(ras_fig3_alcyon_poly$alcyon_ras_pred)) %>%
  # Add alcyon presence points
  addCircles(data = gbif_alcyon_filt_sf, color = "black", weight = 1, radius = 300,
             fillColor = "royalblue", fillOpacity = 0.8,
             group = "A. digitatum presence points"
  ) %>% 
  # Add seafan raster layer
  # addRasterImage(ras_fig3$seafan_ras_pred, colors = col_bin, method = "ngb",
  #                opacity = 0.80, group = "Eunicella verrucosa") %>% 
  # Add alcyon raster layer
  # addRasterImage(ras_fig3$alcyon_ras_pred, colors = col_bin, method = "ngb",
  #                opacity = 0.80, group = "Alcyonium digitatum") %>% 
  # Add layers control
  addLayersControl(options = layersControlOptions(collapsed = FALSE),
                   baseGroups = c("Eunicella verrucosa","Alcyonium digitatum"),
                   overlayGroups = c("E. verrucosa presence points", "A. digitatum presence points")
  ) %>%   
  hideGroup(c("Alcyonium digitatum", "E. verrucosa presence points", "A. digitatum presence points")) %>%
  # Add legend
  addLegend(title = "Habitat suitability", pal = col_bin, values = ras_bin, opacity = 0.80)
l1
saveWidget(l1, file = "figures/Figure3_interactive.html", selfcontained = TRUE)


#--------------#
#
# Figure S3 ####
#
#--------------#

# Plot all probabilities for habitat suitability
ras_figS3 = stack(seafan_pred$seafan_ras_pred,
                  alcyon_pred$alcyon_ras_pred)

# Figure S3
figS3 = tm_shape(ras_figS3)+ tm_raster(title = "Habitat suitability", palette = "-RdYlBu")+
  tm_shape(uk, bbox = bb)+ tm_polygons(lwd = 0.2, border.col = "black")+
  tm_compass(type = "arrow", position = c(0.01,0.08), size = 0.7)+
  tm_scale_bar(position = c(0.01,0.01), breaks = c(0,100,200), lwd = 0.5)+
  tm_layout(panel.labels = fig3_facet_lab, panel.label.size = 0.8,
            legend.outside = FALSE, legend.height = 0.15, legend.text.size = 0.7, legend.title.size = 1)

figS3
tmap_save(tm = figS3, filename = "figures/FigureS3.png", dpi = 600)
tmap_save(tm = figS3, filename = "figures/FigureS3.pdf")



#--------------#
#
# Figure 4 ####
#
#--------------#

# Import raster predictions
load("data/seafan_pred_ras.RData")
seafan_pred
load("data/alcyon_pred_ras.RData")
alcyon_pred

# Function that returns raster predictions showing only probabilities > 0.5
ras_pred_filt = function(raster_file){
  tmpfilter = raster_file < 0.5
  return(mask(raster_file, tmpfilter, maskvalue = 1))
}

# Create raster stack
ras_fig4 = stack(ras_pred_filt(seafan_pred$seafan_ras_pred),
                 ras_pred_filt(seafan_pred$seafan_pred_2050_RCP45),
                 ras_pred_filt(seafan_pred$seafan_pred_2100_RCP45)
)

# Facet labels
fig4_facet_lab = c(expression(italic("Eunicella verrucosa")*":"~"present-day"),"2050 RCP 4.5","2100 RCP 4.5")

# Figure 4 (RCP 4.5)
fig4 = tm_shape(ras_fig4)+ tm_raster(title = "Habitat suitability", palette = "-RdYlBu")+ tm_facets(nrow = 1)+
  tm_shape(uk, bbox = bb)+ tm_polygons(lwd = 0.2, border.col = "black")+
  tm_compass(type = "arrow", position = c(0.01,0.08), size = 0.7)+
  tm_scale_bar(position = c(0.01,0.01), breaks = c(0,100,200), lwd = 0.5)+
  tm_layout(panel.labels = fig4_facet_lab, panel.label.size = 0.8,
            legend.outside = FALSE, legend.height = 0.15, legend.text.size = 0.7, legend.title.size = 1)
fig4
tmap_save(tm = fig4, filename = "figures/Figure4.png", dpi = 600)
tmap_save(tm = fig4, filename = "figures/Figure4.pdf")

# Find difference between predicted and present-day values and plot maps
fig4_diff = stack(ras_fig4$seafan_pred_2050_RCP45 - ras_fig4$seafan_ras_pred,
                  ras_fig4$seafan_pred_2100_RCP45 - ras_fig4$seafan_ras_pred)
plot(fig4_diff)


#--------------#
#
# Figure 4: Interactive Version ####
#
#--------------#

# Function that returns raster predictions showing only probabilities > 0.5
ras_pred_filt = function(raster_file){
  tmpfilter = raster_file < 0.5
  return(mask(raster_file, tmpfilter, maskvalue = 1))
}

# Create raster stack
ras_fig4 = stack(ras_pred_filt(seafan_pred$seafan_ras_pred),
                 ras_pred_filt(seafan_pred$seafan_pred_2050_RCP45),
                 ras_pred_filt(seafan_pred$seafan_pred_2100_RCP45)
                 )

# Convert rasters to polygons
ras_fig4_seafan_pres = rasterToPolygons(ras_fig4$seafan_ras_pred)
ras_fig4_seafan_2050rcp45 = rasterToPolygons(ras_fig4$seafan_pred_2050_RCP45)
ras_fig4_seafan_2100rcp45 = rasterToPolygons(ras_fig4$seafan_pred_2100_RCP45)

# Round habitat suitability to 3 decimal places
round_probs = function(x) { sprintf("%.3f", x) }
ras_fig4_seafan_pres$seafan_ras_pred = round_probs(ras_fig4_seafan_pres$seafan_ras_pred)
ras_fig4_seafan_2050rcp45$seafan_pred_2050_RCP45 = round_probs(ras_fig4_seafan_2050rcp45$seafan_pred_2050_RCP45)
ras_fig4_seafan_2100rcp45$seafan_pred_2100_RCP45 = round_probs(ras_fig4_seafan_2100rcp45$seafan_pred_2100_RCP45)

# Figure 3 interactive version
l2 = leaflet() %>%
  # Boundary
  # fitBounds(lng1 = bb@xmin, lat1 = bb@ymin, lng2 = bb@xmax, lat2 = bb@ymax) %>%
  # Add mouse coordinates at top of map
  addMouseCoordinates() %>%
  # Add an inset minimap
  addMiniMap(
    position = "topright",
    tiles = providers$OpenStreetMap,
    toggleDisplay = TRUE) %>%
  # Add logo and link at the top of the map
  leafem::addLogo(
    img = "https://tomjenkins.netlify.app/project/pinkseafan/featured_hu6b938d805fc03dfee7e876a850589a8a_2399056_720x0_resize_q90_lanczos.JPG",
    url = "https://tomjenkins.netlify.app/project/pinkseafan/",
    src = "remote",
    position = "topleft", width = 150, height = 105
  ) %>% 
  # Add a scalebar
  addScaleBar(
    position = "bottomright",
    options = scaleBarOptions(imperial = FALSE)
  ) %>%
  # Reset map to default setting
  addResetMapButton() %>% 
  # Add measurement tool
  addMeasure(position = "topleft",
             primaryLengthUnit = "meters",
             secondaryLengthUnit = "kilometers",
             primaryAreaUnit = "sqmeters") %>% 
  # Add basemap
  addProviderTiles(providers$OpenStreetMap) %>%
  # Add seafan polygons
  addPolygons(data = ras_fig4_seafan_pres, group = "E. verrucosa present-day",
              fillColor = ~col_bin(ras_fig4_seafan_pres$seafan_ras_pred %>% as.numeric),
              fillOpacity = 0.80,
              stroke = 1, weight = 0.2, color = "black",
              popup = ~htmlEscape(ras_fig4_seafan_pres$seafan_ras_pred)) %>%
  # Add seafan polygons 2050
  addPolygons(data = ras_fig4_seafan_2050rcp45, group = "E. verrucosa 2050",
              fillColor = ~col_bin(ras_fig4_seafan_2050rcp45$seafan_pred_2050_RCP45 %>% as.numeric),
              fillOpacity = 0.80,
              stroke = 1, weight = 0.2, color = "black",
              popup = ~htmlEscape(ras_fig4_seafan_2050rcp45$seafan_pred_2050_RCP45)) %>%
  # Add seafan polygons 2100
  addPolygons(data = ras_fig4_seafan_2100rcp45, group = "E. verrucosa 2100",
              fillColor = ~col_bin(ras_fig4_seafan_2100rcp45$seafan_pred_2100_RCP45 %>% as.numeric),
              fillOpacity = 0.80,
              stroke = 1, weight = 0.2, color = "black",
              popup = ~htmlEscape(ras_fig4_seafan_2100rcp45$seafan_pred_2100_RCP45)) %>%
  # Add seafan presence points
  addCircles(data = gbif_seafan_filt_sf, color = "black", weight = 1, radius = 300,
             fillColor = "deeppink", fillOpacity = 0.8,
             group = "E. verrucosa presence points"
  ) %>% 
  # Add layers control
  addLayersControl(options = layersControlOptions(collapsed = FALSE),
                   baseGroups = c("E. verrucosa present-day",
                                  "E. verrucosa 2050",
                                  "E. verrucosa 2100"),
                   overlayGroups = c("E. verrucosa presence points")
  ) %>%   
  hideGroup(c("E. verrucosa presence points", "E. verrucosa 2050", "E. verrucosa 2100")) %>%
  # Add legend
  addLegend(title = "Habitat suitability", pal = col_bin, values = ras_bin, opacity = 0.80)
l2
saveWidget(l2, file = "figures/Figure4_interactive.html", selfcontained = TRUE)



#--------------#
#
# Figure S4 ####
#
#--------------#

# Figure S4 (RCP 8.5)
ras_figS4 = stack(ras_pred_filt(seafan_pred$seafan_ras_pred),
                  ras_pred_filt(seafan_pred$seafan_pred_2050_RCP85),
                  ras_pred_filt(seafan_pred$seafan_pred_2100_RCP85)
)
figS4_facet_lab = c(expression(italic("Eunicella verrucosa")*":"~"present-day"),"2050 RCP 8.5","2100 RCP 8.5")
figS4 = tm_shape(ras_figS4)+ tm_raster(title = "Habitat suitability", palette = "-RdYlBu")+ tm_facets(nrow = 1)+
  tm_shape(uk, bbox = bb)+ tm_polygons(lwd = 0.2, border.col = "black")+
  tm_compass(type = "arrow", position = c(0.01,0.08), size = 0.7)+
  tm_scale_bar(position = c(0.01,0.01), breaks = c(0,100,200), lwd = 0.5)+
  tm_layout(panel.labels = figS4_facet_lab, panel.label.size = 0.8,
            legend.outside = FALSE, legend.height = 0.15, legend.text.size = 0.7, legend.title.size = 1)
figS4
tmap_save(tm = figS4, filename = "figures/FigureS4.png", dpi = 600)
tmap_save(tm = figS4, filename = "figures/FigureS4.pdf")



#--------------#
#
# Figure S1 ####
#
#--------------#

# Raster predictors to plot
ras_to_plot = c("BO21_tempmin_bdmean","MS_bathy_5m","MS_biogeo06_bathy_slope_5m","BO_calcite",
                "BO_ph","Rock50cm","OrbitalVelMean","TidalVelMean")
ras_predictors = raster::subset(rasters, ras_to_plot)

# Plot raster predictor heatmaps
figS1 = tm_shape(ras_predictors)+
  tm_raster(title = "")+
  tm_facets(free.scales = TRUE)+
  tm_layout(legend.position = c("right","bottom"), legend.height = 0.3, legend.width = 0.3,
            panel.label.size = 0.9)
figS1
tmap_save(tm = figS1, filename = "figures/FigureS1.png", width = 7, height = 4.5, dpi = 900)
tmap_save(tm = figS1, filename = "figures/FigureS1.pdf", width = 7, height = 4.5)


#--------------#
#
# Figure S2 ####
#
#--------------#

# Import background point data
seafan_bgr_pts = read.csv("data/seafan_bgr_pts.csv")
alcyon_bgr_pts = read.csv("data/alcyon_bgr_pts.csv")

# Plot background points
psf_bgpts = ggplot()+
  geom_sf(data = uk, colour = "black", fill = "grey", size = 0.2)+
  geom_point(data = seafan_bgr_pts, aes(x = lon, y = lat),
             shape = 21, fill = "white", colour = "black", size = 1)+
  coord_sf(xlim = c(bb@xmin,bb@xmax), ylim = c(bb@ymin,bb@ymax), expand = FALSE)+
  xlab("Longitude")+
  ylab("Latitude")+
  ggtitle(expression(italic("Eunicella verrucosa")*": background points"))
psf_bgpts
dmf_bgpts = ggplot()+
  geom_sf(data = uk, colour = "black", fill = "grey", size = 0.2)+
  geom_point(data = alcyon_bgr_pts, aes(x = lon, y = lat),
             shape = 21, fill = "white", colour = "black", size = 1)+
  coord_sf(xlim = c(bb@xmin,bb@xmax), ylim = c(bb@ymin,bb@ymax), expand = FALSE)+
  xlab("Longitude")+
  ylab("")+
  ggtitle(expression(italic("Alcyonium digitatum")*": background points"))
dmf_bgpts

# Figure S2
figS2 = ggarrange(psf_bgpts + fig1_theme,
                  dmf_bgpts + fig1_theme)
figS2
ggsave(plot = figS2, file = "figures/FigureS2.png", width = 8, height = 5, dpi = 600)
ggsave(plot = figS2, file = "figures/FigureS2.pdf", width = 8, height = 5)



