# --------------------------- #
#
# Species Distribution Modelling
#
# Jenkins and Stevens (2021) in prep
#
# R script purpose:
# 1. Pink sea fan: Maxent SDM
# 2. Dead man's fingers: Maxent SDM
#
# Author: Tom Jenkins
# Email: tom.l.jenkins@outlook.com
# Website: https://tomjenkins.netlify.app
#
# --------------------------- #

# Load libraries
# .libPaths("C:/R-4.1.0/library/")
library(tidyverse)
library(raster)
library(tmap)
library(tmaptools)
library(sf)
library(ENMeval)
library(maxnet)
library(rJava)
library(dismo)
library(usdm)


# --------------- #
#
# Collinearity of predictor variables ####
#
# --------------- #

# Import raster predictor variables
load("data/raster_predictors.RData")

# Subset present-day variables from raster stack by removing all variable names with "RCP"
ras_predictors = raster::subset(rasters, grep("RCP", names(rasters), value = TRUE, invert = TRUE))
names(ras_predictors)

# Assess collinearity of raster predictor variables
ras_vif = vifcor(ras_predictors, maxobservations = ncell(ras_predictors), th = 1)
ras_vif
ras_cor = ras_vif@corMatrix

# Vector of variables to remove from raster stack
variables_to_remove = c("BO21_tempmin_ss","BO21_ppmean_bdmean")
ras_predictors = dropLayer(ras_predictors, variables_to_remove)
names(ras_predictors)

# Re-check variance inflation factor for rasters
vifcor(ras_predictors, maxobservations = ncell(ras_predictors), th = 1)


# ================= #
#
# 1. Pink sea fan: Maxent SDM ####
#
# ================= #

# Import pink sea fan data
seafan_data = read.csv("data/sdm_seafan_input_data.csv")

# Subset occurrence points
seafan_occur_pts = seafan_data %>%
  dplyr::filter(occur.backgr == 1) %>%
  dplyr::select(Lon, Lat)

# Subset background points
seafan_backgr_pts = seafan_data %>%
  dplyr::filter(occur.backgr == 0) %>%
  dplyr::select(Lon, Lat)

# Visualise the 'block' spatial partition method
seafan_block = get.block(seafan_occur_pts, seafan_backgr_pts, orientation = "lat_lat")
table(seafan_block$occs.grp)
evalplot.grps(pts = seafan_occur_pts, pts.grp = seafan_block$occs.grp, envs = ras_predictors$MS_bathy_5m) + 
  ggplot2::ggtitle("Spatial block partitions: presences")


#--------------#
#
# ENMevaluate Maxent ####
#
#--------------#

# Run MaxEnt using maxent.jar from dismo package
set.seed(123)
mx_seafan = ENMevaluate(occs = seafan_occur_pts,
                        envs = ras_predictors,
                        bg = seafan_backgr_pts,
                        algorithm = "maxent.jar",
                        partitions = "block",
                        partition.settings = list(orientation = "lat_lat"),
                        tune.args = list(fc = c("L","Q","LQ","H"),
                                         rm = 1:5),
                        parallel = TRUE,
                        numCores = 4
                        )
mx_seafan

# ENMevaluate overall results table
mx_seafan@results

# ENMevaluate partitions results table
eval.partition.method(mx_seafan)
eval.results.partitions(mx_seafan)

# Model selection: lowest AICc score
# In practice, models with delta AICc scores less than 2 are usually considered 
# statistically equivalent
mx_seafan@results %>% arrange(delta.AICc) %>% dplyr::select(tune.args, delta.AICc)
opt_score_seafan = mx_seafan@results %>% filter(delta.AICc == 0)
opt_score_seafan

# Model selection: sequential criteria
# First, selects models with lowest average test omission rate
# Second, select models with highest average validation AUC
# opt_score_seafan = mx_seafan@results %>%
#   dplyr::filter(or.10p.avg == min(or.10p.avg)) %>%
#   dplyr::filter(auc.val.avg == max(auc.val.avg))
# opt_score_seafan

# Extract optimum model from ENMevaluation object
opt_model_seafan = eval.models(mx_seafan)[[opt_score_seafan$tune.args]]

# Marginal response curves for optimal model
dismo::response(opt_model_seafan)

# Extract variable importance from ENMevaluation object
seafan_var = mx_seafan@variable.importance[opt_score_seafan$tune.args][[1]] %>%
  arrange(percent.contribution)
seafan_var
write_csv(seafan_var, file = "data/seafan_var.csv")

# Extract raster predictions from ENMevaluation object
seafan_ras_pred = eval.predictions(mx_seafan)[[opt_score_seafan$tune.args]]

# Plot predictions showing probability of presence
tmS1 = tm_shape(seafan_ras_pred)+ tm_raster()

# Function that returns raster predictions showing only probabilities > 0.5
ras_pred_filt = function(raster_file){
  tmpfilter = raster_file < 0.5
  return(mask(raster_file, tmpfilter, maskvalue = 1))
}
tmS2 = tm_shape(ras_pred_filt(seafan_ras_pred))+ tm_raster()
tmap_arrange(tmS1, tmS2)

# Execute Null model
# Run null simulations with 100 iterations to get a reasonable null distribution 
# for comparisons with the empirical values
# set.seed(123)
# seafan_null = ENMnulls(mx_seafan,
#                        mod.settings = list(fc = as.character(opt_score_seafan$fc),
#                                            rm = as.numeric(opt_score_seafan$rm)),
#                        no.iter = 100,
#                        parallel = TRUE,
#                        numCores = 4)
# save(seafan_null, file = "data/seafan_null_model.RData")
load("data/seafan_null_model.RData")

# Comparison between the empirical and simulated results
null.emp.results(seafan_null)

# Plots null model results
evalplot.nulls(seafan_null, stats = c("or.10p", "auc.val"), plot.type = "histogram")
evalplot.nulls(seafan_null, stats = c("or.10p", "auc.val"), plot.type = "violin")


#--------------#
#
# Future predictions ####
#
#--------------#

# Create new raster predictor object without present temperature layers
seafan_2050_RCP45 = names(ras_predictors)[2:nlayers(ras_predictors)]
seafan_2050_RCP45 = raster::subset(ras_predictors, seafan_2050_RCP45)
names(seafan_2050_RCP45)

# Append 2050 RCP 4.5 future temperature layers to raster object
seafan_2050_RCP45 = raster::stack(rasters$BO21_RCP45_2050_tempmin_bdmean,
                                  seafan_2050_RCP45)
names(seafan_2050_RCP45)

# Change 2050 RCP 4.5 future temperature names to match variable names in model
names(seafan_2050_RCP45)[1] = names(ras_predictors)[1]
names(seafan_2050_RCP45)

# Maxent model predictions with 2050 RCP 4.5 future raster data
seafan_pred_2050_RCP45 = raster::predict(model = opt_model_seafan, object = seafan_2050_RCP45)

# Repeat for 2100 RCP 4.5 future temperature layers
seafan_2100_RCP45 = names(ras_predictors)[2:nlayers(ras_predictors)]
seafan_2100_RCP45 = raster::subset(ras_predictors, seafan_2100_RCP45)
names(seafan_2100_RCP45)
seafan_2100_RCP45 = raster::stack(rasters$BO21_RCP45_2100_tempmin_bdmean,
                                  seafan_2100_RCP45)
names(seafan_2100_RCP45)
names(seafan_2100_RCP45)[1] = names(ras_predictors)[1]
names(seafan_2100_RCP45)
seafan_pred_2100_RCP45 = raster::predict(model = opt_model_seafan, object = seafan_2100_RCP45)


# Repeat for 2050 RCP 8.5 future temperature layers
seafan_2050_RCP85 = names(ras_predictors)[2:nlayers(ras_predictors)]
seafan_2050_RCP85 = raster::subset(ras_predictors, seafan_2050_RCP85)
names(seafan_2050_RCP85)
seafan_2050_RCP85 = raster::stack(rasters$BO21_RCP85_2050_tempmin_bdmean,
                                  seafan_2050_RCP85)
names(seafan_2050_RCP85)
names(seafan_2050_RCP85)[1] = names(ras_predictors)[1]
names(seafan_2050_RCP85)
seafan_pred_2050_RCP85 = raster::predict(model = opt_model_seafan, object = seafan_2050_RCP85)

# Repeat for 2100 RCP 8.5 future temperature layers
seafan_2100_RCP85 = names(ras_predictors)[2:nlayers(ras_predictors)]
seafan_2100_RCP85 = raster::subset(ras_predictors, seafan_2100_RCP85)
names(seafan_2100_RCP85)
seafan_2100_RCP85 = raster::stack(rasters$BO21_RCP85_2100_tempmin_bdmean,
                                  seafan_2100_RCP85)
names(seafan_2100_RCP85)
names(seafan_2100_RCP85)[1] = names(ras_predictors)[1]
names(seafan_2100_RCP85)
seafan_pred_2100_RCP85 = raster::predict(model = opt_model_seafan, object = seafan_2100_RCP85)


#--------------#
#
# Visualise future predictions ####
#
#--------------#

# Present-day
tmap_arrange(tmS1, tmS2)

# RCP 4.5
# tmS3 = tm_shape(seafan_pred_2050_RCP45)+ tm_raster()
# tmS4 = tm_shape(seafan_pred_2100_RCP45)+ tm_raster()
tm_shape(ras_pred_filt(seafan_pred_2050_RCP45))+ tm_raster()
tm_shape(ras_pred_filt(seafan_pred_2100_RCP45))+ tm_raster()

# RCP 8.5
# tmS5 = tm_shape(seafan_pred_2050_RCP85)+ tm_raster()
# tmS6 = tm_shape(seafan_pred_2100_RCP85)+ tm_raster()
tm_shape(ras_pred_filt(seafan_pred_2050_RCP85))+ tm_raster()
tm_shape(ras_pred_filt(seafan_pred_2100_RCP85))+ tm_raster()

# Export raster predictions
seafan_pred = raster::stack(seafan_ras_pred, seafan_pred_2050_RCP45, seafan_pred_2100_RCP45, seafan_pred_2050_RCP85, seafan_pred_2100_RCP85)
names(seafan_pred) = c("seafan_ras_pred", "seafan_pred_2050_RCP45", "seafan_pred_2100_RCP45", "seafan_pred_2050_RCP85", "seafan_pred_2100_RCP85")
save(seafan_pred, file = "data/seafan_pred_ras.RData")


# ================= #
#
# 2. Dead man's fingers: Maxent SDM ####
#
# ================= #

# Import dead man's fingers data
alcyon_data = read.csv("data/sdm_alcyon_input_data.csv")

# Subset occurrence points
alcyon_occur_pts = alcyon_data %>%
  dplyr::filter(occur.backgr == 1) %>%
  dplyr::select(Lon, Lat)

# Subset background points
alcyon_backgr_pts = alcyon_data %>%
  dplyr::filter(occur.backgr == 0) %>%
  dplyr::select(Lon, Lat)

# Visualise the 'block' spatial partition method
alcyon_block = get.block(alcyon_occur_pts, alcyon_backgr_pts, orientation = "lat_lat")
table(alcyon_block$occs.grp)
evalplot.grps(pts = alcyon_occur_pts, pts.grp = alcyon_block$occs.grp, envs = ras_predictors$MS_bathy_5m) + 
  ggplot2::ggtitle("Spatial block partitions: presences")


#--------------#
#
# ENMevaluate Maxent ####
#
#--------------#

# Run Maxent using maxent.jar from dismo package
set.seed(123)
mx_alcyon = ENMevaluate(occs = alcyon_occur_pts,
                        envs = ras_predictors,
                        bg = alcyon_backgr_pts,
                        algorithm = "maxent.jar",
                        partitions = "block",
                        partition.settings = list(orientation = "lat_lat"),
                        tune.args = list(fc = c("L","Q","LQ","H"),
                                         rm = 1:5),
                        parallel = TRUE,
                        numCores = 4
)
mx_alcyon

# ENMevaluate overall results table
mx_alcyon@results

# ENMevaluate partitions results table
eval.partition.method(mx_alcyon)
eval.results.partitions(mx_alcyon)

# Model selection: lowest AICc score
# In practice, models with delta AICc scores less than 2 are usually considered 
# statistically equivalent
mx_alcyon@results %>% arrange(delta.AICc) %>% dplyr::select(tune.args, delta.AICc)
opt_score_alcyon = mx_alcyon@results %>% filter(delta.AICc == 0)
opt_score_alcyon

# Model selection: sequential criteria
# First, selects models with lowest average test omission rate
# Second, select models with highest average validation AUC
# opt_score_alcyon = mx_alcyon@results %>%
#   dplyr::filter(or.10p.avg == min(or.10p.avg)) %>%
#   dplyr::filter(auc.val.avg == max(auc.val.avg))
# opt_score_alcyon

# Extract optimal model from ENMevaluation object
opt_model_alcyon = eval.models(mx_alcyon)[[opt_score_alcyon$tune.args]]

# Marginal response curves for optimal model
dismo::response(opt_model_alcyon)

# Extract variable importance from ENMevaluation object
alcyon_var = mx_alcyon@variable.importance[opt_score_alcyon$tune.args][[1]] %>%
  arrange(percent.contribution)
alcyon_var
write_csv(alcyon_var, file = "data/alcyon_var.csv")

# Extract raster predictions from ENMevaluation object
alcyon_ras_pred = eval.predictions(mx_alcyon)[[opt_score_alcyon$tune.args]]

# Plot predictions showing probability of presence
tmA1 = tm_shape(alcyon_ras_pred)+ tm_raster()

# Plot predictions showing probabilities > 0.5
tmA2 = tm_shape(ras_pred_filt(alcyon_ras_pred))+ tm_raster()
tmap_arrange(tmA1, tmA2)

# Execute Null model
# Run null simulations with 100 iterations to get a reasonable null distribution 
# for comparisons with the empirical values
# set.seed(123)
# alcyon_null = ENMnulls(mx_alcyon,
#                        mod.settings = list(fc = as.character(opt_score_alcyon$fc),
#                                            rm = as.numeric(opt_score_alcyon$rm)),
#                        no.iter = 100,
#                        parallel = TRUE,
#                        numCores = 4)
# save(alcyon_null, file = "data/alcyon_null_model.RData")
load("data/alcyon_null_model.RData")

# Comparison between the empirical and simulated results
null.emp.results(alcyon_null)

# Plots null model results
evalplot.nulls(alcyon_null, stats = c("or.10p", "auc.val"), plot.type = "histogram")
evalplot.nulls(alcyon_null, stats = c("or.10p", "auc.val"), plot.type = "violin")


#--------------#
#
# Future predictions ####
#
#--------------#

# Create new raster predictor object without present temperature layers
alcyon_2050_RCP45 = names(ras_predictors)[2:nlayers(ras_predictors)]
alcyon_2050_RCP45 = raster::subset(ras_predictors, alcyon_2050_RCP45)
names(alcyon_2050_RCP45)

# Append 2050 RCP 4.5 future temperature layers to raster object
alcyon_2050_RCP45 = raster::stack(rasters$BO21_RCP45_2050_tempmin_bdmean,
                                  alcyon_2050_RCP45)
names(alcyon_2050_RCP45)

# Change 2050 RCP 4.5 future temperature names to match variable names in model
names(alcyon_2050_RCP45)[1] = names(ras_predictors)[1]
names(alcyon_2050_RCP45)

# Maxent model predictions with 2050 RCP 4.5 future raster data
alcyon_pred_2050_RCP45 = raster::predict(model = opt_model_alcyon, object = alcyon_2050_RCP45)

# Repeat for 2100 RCP 4.5 future temperature layers
alcyon_2100_RCP45 = names(ras_predictors)[2:nlayers(ras_predictors)]
alcyon_2100_RCP45 = raster::subset(ras_predictors, alcyon_2100_RCP45)
names(alcyon_2100_RCP45)
alcyon_2100_RCP45 = raster::stack(rasters$BO21_RCP45_2100_tempmin_bdmean,
                                  alcyon_2100_RCP45)
names(alcyon_2100_RCP45)
names(alcyon_2100_RCP45)[1] = names(ras_predictors)[1]
names(alcyon_2100_RCP45)
alcyon_pred_2100_RCP45 = raster::predict(model = opt_model_alcyon, object = alcyon_2100_RCP45)


# Repeat for 2050 RCP 8.5 future temperature layers
alcyon_2050_RCP85 = names(ras_predictors)[2:nlayers(ras_predictors)]
alcyon_2050_RCP85 = raster::subset(ras_predictors, alcyon_2050_RCP85)
names(alcyon_2050_RCP85)
alcyon_2050_RCP85 = raster::stack(rasters$BO21_RCP85_2050_tempmin_bdmean,
                                  alcyon_2050_RCP85)
names(alcyon_2050_RCP85)
names(alcyon_2050_RCP85)[1] = names(ras_predictors)[1]
names(alcyon_2050_RCP85)
alcyon_pred_2050_RCP85 = raster::predict(model = opt_model_alcyon, object = alcyon_2050_RCP85)

# Repeat for 2100 RCP 8.5 future temperature layers
alcyon_2100_RCP85 = names(ras_predictors)[2:nlayers(ras_predictors)]
alcyon_2100_RCP85 = raster::subset(ras_predictors, alcyon_2100_RCP85)
names(alcyon_2100_RCP85)
alcyon_2100_RCP85 = raster::stack(rasters$BO21_RCP85_2100_tempmin_bdmean,
                                  alcyon_2100_RCP85)
names(alcyon_2100_RCP85)
names(alcyon_2100_RCP85)[1] = names(ras_predictors)[1]
names(alcyon_2100_RCP85)
alcyon_pred_2100_RCP85 = raster::predict(model = opt_model_alcyon, object = alcyon_2100_RCP85)


#--------------#
#
# Visualise future predictions ####
#
#--------------#

# Present-day
tmap_arrange(tmA1, tmA2)

# RCP 4.5
# tmA3 = tm_shape(alcyon_pred_2050_RCP45)+ tm_raster()
# tmA4 = tm_shape(alcyon_pred_2100_RCP45)+ tm_raster()
tm_shape(ras_pred_filt(alcyon_pred_2050_RCP45))+ tm_raster()
tm_shape(ras_pred_filt(alcyon_pred_2100_RCP45))+ tm_raster()

# RCP 8.5
# tmA5 = tm_shape(alcyon_pred_2050_RCP85)+ tm_raster()
# tmA6 = tm_shape(alcyon_pred_2100_RCP85)+ tm_raster()
tm_shape(ras_pred_filt(alcyon_pred_2050_RCP85))+ tm_raster()
tm_shape(ras_pred_filt(alcyon_pred_2100_RCP85))+ tm_raster()

# Export raster predictions
alcyon_pred = raster::stack(alcyon_ras_pred, alcyon_pred_2050_RCP45, alcyon_pred_2100_RCP45, alcyon_pred_2050_RCP85, alcyon_pred_2100_RCP85)
names(alcyon_pred) = c("alcyon_ras_pred", "alcyon_pred_2050_RCP45", "alcyon_pred_2100_RCP45", "alcyon_pred_2050_RCP85", "alcyon_pred_2100_RCP85")
save(alcyon_pred, file = "data/alcyon_pred_ras.RData")

