# --------------------------- #
#
# Species Distribution Modelling
#
# Jenkins and Stevens (2022) in prep
#
# R script purpose:
# 1. Pink sea fan: Maxent SDM
# 2. Dead man's fingers: Maxent SDM
#
# Author: Tom Jenkins
# Email: tom.l.jenkins@outlook.com
#
# --------------------------- #

# Load libraries
library(tidyverse)
library(raster)
library(terra)
library(tmap)
library(tmaptools)
library(sf)
library(ENMeval)
library(maxnet)
library(rJava)
library(dismo)
library(usdm)

# Import raster predictor variables
bathymetry = raster("../data/raster_predictors/bathymetry.tif")
slope = raster("../data/raster_predictors/slope.tif")
strathclyde = raster::stack("../data/raster_predictors/strathclyde.tif")
env_ras = raster::stack("../data/raster_predictors/env_rasters.tif")
# biooracle = raster::stack("../data/raster_predictors/biooracle.tif")

# Subset present-day variables from raster stack
ras_predictors = raster::stack(bathymetry, slope, strathclyde, env_ras) %>% 
  raster::subset(., str_subset(names(.), "2081_2100", negate = TRUE))
names(ras_predictors)

# Assess collinearity of raster predictor variables
ras_vif = vifcor(ras_predictors, maxobservations = ncell(ras_predictors), th = 1)
ras_vif
ras_cor = ras_vif@corMatrix

# Vector of variables to remove from raster stack
variables_to_remove = c("bathymetry",
                        "ph_FromKrige_3km_1951_2000",
                        "Arag_FromKrige_3km_1951_2000",
                        "Epc_FromKrige_3km_1951_2000")
ras_predictors = dropLayer(ras_predictors, variables_to_remove)
names(ras_predictors)

# Re-check variance inflation factor for rasters
vifcor(ras_predictors, maxobservations = ncell(ras_predictors), th = 0.70)

# Export raster predictors as .RData object
raster::writeRaster(ras_predictors, filename = "../data/ras_predictors.grd")


# --------------- #
#
# 1. Pink sea fan: Maxent SDM ####
#
# --------------- #

# Import pink sea fan data
psf_data = read_csv("../data/sdm_input_data.csv") %>% 
  dplyr::filter(species == "E_verrucosa")
psf_data

# Assess collinearity of raster predictor variables
# env_vars_to_remove = c()
# env_vars_to_remove = variables_to_remove
# psf_data %>%
#   dplyr::select(-c(1:5, str_which(names(psf_data), "2100"), env_vars_to_remove)) %>%
#   as.data.frame %>%
#   usdm::vifcor(., th = 1)

# Remove correlated variables from raster stack
# ras_predictors = dropLayer(ras_predictors, env_vars_to_remove)
# names(ras_predictors)

# Subset occurrence points
psf_occur_pts = psf_data %>%
  dplyr::filter(pres_backgr == 1) %>%
  dplyr::select(lon, lat)

# Subset background points
psf_backgr_pts = psf_data %>%
  dplyr::filter(pres_backgr == 0) %>%
  dplyr::select(lon, lat)


# --------------- #
#
# Block spatial partition (lat-lat) ####
#
# --------------- #

# Visualise the 'block' spatial partition method
psf_block = get.block(psf_occur_pts, psf_backgr_pts, orientation = "lat_lat")
table(psf_block$occs.grp)
evalplot.grps(pts = psf_occur_pts, pts.grp = psf_block$occs.grp, envs = bathymetry) + 
  ggplot2::ggtitle("Spatial block partitions: presences")
evalplot.grps(pts = psf_backgr_pts, pts.grp = psf_block$bg.grp, envs = bathymetry) + 
  ggplot2::ggtitle("Spatial block partitions: background")


# --------------- #
#
# ENMevaluate Maxent ####
#
# --------------- #

# Run MaxEnt using maxent.jar from dismo package
set.seed(123)
psf_maxent = ENMevaluate(
  occs = psf_occur_pts,
  envs = ras_predictors,
  bg = psf_backgr_pts,
  algorithm = "maxent.jar",
  partitions = "block",
  partition.settings = list(orientation = "lat_lat"),
  tune.args = list(fc = c("L","LQ","H"), rm = 1:5),
  parallel = TRUE,
  numCores = 4
)
psf_maxent

# ENMevaluate overall results table
# psf_maxent@results

# ENMevaluate partitions results table
# eval.partition.method(psf_maxent)
# eval.results.partitions(psf_maxent)
# psf_maxent@results %>% arrange(fc, rm) %>% dplyr::select(tune.args, delta.AICc)
psf_maxent@results %>% dplyr::select(fc, rm, auc.val.avg, auc.val.sd, cbi.val.avg, cbi.val.sd)

# Model selection: lowest AICc score
# In practice, models with delta AICc scores less than 2 are usually considered 
# statistically equivalent
psf_maxent@results %>% arrange(delta.AICc) %>% dplyr::select(tune.args, delta.AICc)
psf_opt_score = psf_maxent@results %>%
  arrange(delta.AICc) %>% 
  slice(1)
psf_opt_score

# Model selection: sequential criteria
# First, selects models with lowest average test omission rate
# Second, select models with highest average validation AUC
# psf_opt_score = psf_maxent@results %>%
#   dplyr::filter(or.10p.avg == min(or.10p.avg)) %>%
#   dplyr::filter(auc.val.avg == max(auc.val.avg))
#   dplyr::filter(cbi.val.avg == max(cbi.val.avg))
# psf_opt_score

# Extract optimum model from ENMevaluation object
psf_opt_model = eval.models(psf_maxent)[[psf_opt_score$tune.args]]

# Marginal response curves for optimal model
dismo::response(psf_opt_model)

# Extract variable importance from ENMevaluation object
psf_var = psf_maxent@variable.importance[psf_opt_score$tune.args][[1]] %>%
  arrange(percent.contribution)
psf_var
write_csv(psf_var, file = "../data/pinkseafan_varcontrib.csv")

# Extract predictions as a raster object
psf_pred = eval.predictions(psf_maxent)[[psf_opt_score$tune.args]]

# Plot predictions showing habitat suitability predictions
psf_tm1 = tm_shape(psf_pred)+ tm_raster()
psf_tm1

# Execute Null model
# Run null simulations with 100 iterations to get a reasonable null distribution 
# for comparisons with the empirical values
# set.seed(123)
# psf_null = ENMnulls(
#   psf_maxent,
#   mod.settings = list(
#     fc = as.character(psf_opt_score$fc), rm = as.numeric(psf_opt_score$rm)),
#   no.iter = 100,
#   parallel = TRUE,
#   numCores = 4)
# save(psf_null, file = "../data/pinkseafan_null_model.RData")
load("../data/pinkseafan_null_model.RData")

# Comparison between the empirical and simulated results
null.emp.results(psf_null)

# Plots null model results
evalplot.nulls(psf_null, stats = c("or.10p", "auc.val"), plot.type = "histogram")
# evalplot.nulls(psf_null, stats = c("or.10p", "auc.val"), plot.type = "violin")


# --------------- #
#
# Future predictions ####
#
# --------------- #

# Create new raster predictor object with future raster layers
future_layers = str_subset(names(ras_predictors), "2000") %>% str_replace_all(., "1951_2000", "2081_2100")
ras_future = raster::stack(
  raster::subset(ras_predictors, str_subset(names(ras_predictors), "1951", negate = TRUE)),
  raster::subset(env_ras, future_layers))
names(ras_future)

# Change names of future raster layers to match variable names in model
names(ras_future) = names(ras_predictors)
names(ras_future)

# Maxent model predictions with future raster data
psf_pred_future = raster::predict(
  model = psf_opt_model,
  object = ras_future
)

# Future habitat suitability index
psf_tm2 = tm_shape(psf_pred_future)+ tm_raster()

# Present versus future
tmap_arrange(psf_tm1, psf_tm2)

# Export raster predictions
psf_raster = raster::stack(psf_pred, psf_pred_future) %>% 
  set_names(c("psf_pred","psf_pred_future"))
raster::writeRaster(psf_raster, filename = "../data/pinkseafan_ras.grd")



# --------------- #
#
# 2. Dead man's fingers: Maxent SDM ####
#
# --------------- #

# Import dead man's fingers data
dmf_data = read_csv("../data/sdm_input_data.csv") %>% 
  dplyr::filter(species == "A_digitatum")
dmf_data

# Assess collinearity of raster predictor variables
# env_vars_to_remove = c()
# env_vars_to_remove = variables_to_remove
# dmf_data %>%
#   dplyr::select(-c(1:5, str_which(names(dmf_data), "2100"), env_vars_to_remove)) %>%
#   as.data.frame %>%
#   usdm::vifcor(., th = 1)

# Remove correlated variables from raster stack
# ras_predictors = dropLayer(ras_predictors, env_vars_to_remove)
# names(ras_predictors)

# Subset occurrence points
dmf_occur_pts = dmf_data %>%
  dplyr::filter(pres_backgr == 1) %>%
  dplyr::select(lon, lat)

# Subset background points
dmf_backgr_pts = dmf_data %>%
  dplyr::filter(pres_backgr == 0) %>%
  dplyr::select(lon, lat)


# --------------- #
#
# Block spatial partition (lat-lat) ####
#
# --------------- #

# Visualise the 'block' spatial partition method
dmf_block = get.block(dmf_occur_pts, dmf_backgr_pts, orientation = "lat_lat")
table(dmf_block$occs.grp)
evalplot.grps(pts = dmf_occur_pts, pts.grp = dmf_block$occs.grp, envs = bathymetry) + 
  ggplot2::ggtitle("Spatial block partitions: presences")
evalplot.grps(pts = dmf_backgr_pts, pts.grp = dmf_block$bg.grp, envs = bathymetry) + 
  ggplot2::ggtitle("Spatial block partitions: background")


# --------------- #
#
# ENMevaluate Maxent ####
#
# --------------- #

# Run MaxEnt using maxent.jar from dismo package
set.seed(123)
dmf_maxent = ENMevaluate(
  occs = dmf_occur_pts,
  envs = ras_predictors,
  bg = dmf_backgr_pts,
  algorithm = "maxent.jar",
  partitions = "block",
  partition.settings = list(orientation = "lat_lat"),
  tune.args = list(fc = c("L","LQ","H"), rm = 1:5),
  parallel = TRUE,
  numCores = 4
)
dmf_maxent

# ENMevaluate overall results table
# dmf_maxent@results

# ENMevaluate partitions results table
# eval.partition.method(dmf_maxent)
# eval.results.partitions(dmf_maxent)

# Model selection: lowest AICc score
# In practice, models with delta AICc scores less than 2 are usually considered 
# statistically equivalent
dmf_maxent@results %>% arrange(delta.AICc) %>% dplyr::select(tune.args, delta.AICc) %>% head
dmf_opt_score = dmf_maxent@results %>%
  arrange(delta.AICc) %>% 
  slice(1)
dmf_opt_score

# Model selection: sequential criteria
# First, selects models with lowest average test omission rate
# Second, select models with highest average validation AUC
# dmf_opt_score = dmf_maxent@results %>%
#   dplyr::filter(or.10p.avg == min(or.10p.avg)) %>%
#   # dplyr::filter(auc.val.avg == max(auc.val.avg))
#   dplyr::filter(cbi.val.avg == max(cbi.val.avg))
# dmf_opt_score

# Extract optimum model from ENMevaluation object
dmf_opt_model = eval.models(dmf_maxent)[[dmf_opt_score$tune.args]]

# Marginal response curves for optimal model
dismo::response(dmf_opt_model)

# Extract variable importance from ENMevaluation object
dmf_var = dmf_maxent@variable.importance[dmf_opt_score$tune.args][[1]] %>%
  arrange(percent.contribution)
dmf_var
write_csv(dmf_var, file = "../data/deadseafan_varcontrib.csv")

# Extract predictions as a raster object
dmf_pred = eval.predictions(dmf_maxent)[[dmf_opt_score$tune.args]]

# Plot predictions showing habitat suitability predictions
dmf_tm1 = tm_shape(dmf_pred)+ tm_raster()
dmf_tm1

# Execute Null model
# Run null simulations with 100 iterations to get a reasonable null distribution 
# for comparisons with the empirical values
# set.seed(123)
# dmf_null = ENMnulls(
#   dmf_maxent,
#   mod.settings = list(
#     fc = as.character(dmf_opt_score$fc), rm = as.numeric(dmf_opt_score$rm)),
#   no.iter = 100,
#   parallel = TRUE,
#   numCores = 4)
# save(dmf_null, file = "../data/deadmansfingers_null_model.RData")
load("../data/deadmansfingers_null_model.RData")

# Comparison between the empirical and simulated results
null.emp.results(dmf_null)

# Plots null model results
evalplot.nulls(dmf_null, stats = c("or.10p", "auc.val"), plot.type = "histogram")
evalplot.nulls(dmf_null, stats = c("or.10p", "auc.val"), plot.type = "violin")


# --------------- #
#
# Future predictions ####
#
# --------------- #

# Create new raster predictor object with future raster layers
future_layers = str_subset(names(ras_predictors), "2000") %>% str_replace_all(., "1951_2000", "2081_2100")
ras_future = raster::stack(
  raster::subset(ras_predictors, str_subset(names(ras_predictors), "1951", negate = TRUE)),
  raster::subset(env_ras, future_layers))
names(ras_future)

# Change names of future raster layers to match variable names in model
names(ras_future) = names(ras_predictors)
names(ras_future)

# Maxent model predictions with future raster data
dmf_pred_future = raster::predict(
  model = dmf_opt_model,
  object = ras_future
)

# Future habitat suitability index
dmf_tm2 = tm_shape(dmf_pred_future)+ tm_raster()

# Present versus future
tmap_arrange(dmf_tm1, dmf_tm2)

# Export raster predictions
dmf_raster = raster::stack(dmf_pred, dmf_pred_future) %>% 
  set_names(c("dmf_pred","dmf_pred_future"))
raster::writeRaster(dmf_raster, filename = "../data/deadmansfingers_ras.grd")



# --------------- #
#
# Extract info for figures / supplementary ####
#
# --------------- #

# Table of evaluation stats for all models
library(gt)
plot_gt = function(maxent_res_df, title = ""){
  eval_stats = maxent_res_df %>%
    dplyr::select(fc, rm, auc.val.avg, auc.val.sd, cbi.val.avg, cbi.val.sd, delta.AICc) %>% 
    mutate(fc = factor(fc, levels = c("L","LQ","H"))) %>% 
    arrange(fc, rm) %>% 
    mutate_if(is.numeric, round, 2)
  eval_stats %>%
    gt() %>% 
    tab_header(title = title) %>% 
    tab_style(
      style = list(cell_text(style = "italic")),
      locations = cells_title(groups = "title")) %>% 
    tab_style(
      style = list(cell_text(weight = "bold")),
      locations = cells_body(columns = everything(), rows = 11))
}
(psf_gt = plot_gt(psf_maxent@results, title = "Eunicella verrucosa"))
(dmf_gt = plot_gt(dmf_maxent@results, title = "Alcyonium digitatum"))
gtsave(psf_gt, filename = "../figures/gt_psf.tex")
gtsave(dmf_gt, filename = "../figures/gt_dmf.tex")


