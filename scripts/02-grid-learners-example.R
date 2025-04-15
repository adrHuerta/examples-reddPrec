rm(list = ls())
# https://www.meteoswiss.admin.ch/climate/the-climate-of-switzerland/records-and-extremes/extreme-value-analyses/high-impact-precipitation-events/24-28-july-2014/precipitation-and-temperature.html

library(xts)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(reddPrec)
source("R/metrics.R")

# ch data
ch_data <- readRDS("data/ch_gf_data.RDS")

# ch topo variables
covs_ch_data <- rast("data/topo_data/topo_vars_ch.nc")
names(covs_ch_data) <- c('aspectcosine', 'aspectsine', 'dist2coast', 'dx', 'dxx', 'dy', 'dyy', 'eastness',
                         'elevation', 'latitude', 'longitude', 'northness', 'pcurv', 'roughness', 'slope', 'tcurv', 'tpi',
                         'tri', 'vrm')

# pca of topo variables
pca_cvs_ch_data <- princomp(covs_ch_data[[-c(11, 10, 9)]])
pca_cvs_ch_data <- predict(covs_ch_data[[-c(11, 10, 9)]], pca_cvs_ch_data, index=1:4)
names(pca_cvs_ch_data) <- c("pc1", "pc2", "pc3", "pc4")

#
lon_lat_alt <- covs_ch_data[[c(11, 10, 9)]]
names(lon_lat_alt) <- c("lon", "lat", "alt")
lon_lat_alt <- c(lon_lat_alt, pca_cvs_ch_data) 

# surface reflectance data
modis_sur_refl_b01 <- lapply(dir("data/modis_data", pattern = "MODIS_Aqua_2014",
                                 full.names = TRUE),
                             
                             function(x) {
                               # gap-filling
                               out <- project(rast(x)[[1]], lon_lat_alt)
                               out <- focal(out, w=17, fun=mean, na.policy="only")
                               out <- focal(out, 17, "modal", na.policy="only")
                               out
                               
                             }
)
modis_sur_refl_b01 <- do.call(c, modis_sur_refl_b01)
names(modis_sur_refl_b01) <- rep("sur_refl_b01", 5)

modis_sur_refl_b02 <- lapply(dir("data/modis_data", pattern = "MODIS_Aqua_2014",
                                 full.names = TRUE),
                             function(x) {
                               # gap-filling
                               out <- project(rast(x)[[2]], lon_lat_alt)
                               out <- focal(out, w=17, fun=mean, na.policy="only")
                               out <- focal(out, 17, "modal", na.policy="only")
                               out
                               
                             }
)
modis_sur_refl_b02 <- do.call(c, modis_sur_refl_b02)
names(modis_sur_refl_b02) <- rep("sur_refl_b02", 5)

# setting arguments for gridPcp
dOe <- seq(as.Date("2014-07-24"), as.Date("2014-07-28"), by = "day")

grid <- c(lon_lat_alt)
grid <- crop(grid, ext(6, 10.5, 45.75, 47.9))

dyncovars <- sds(modis_sur_refl_b01, modis_sur_refl_b02)
dyncovars <- crop(dyncovars, ext(6, 10.5, 45.75, 47.9))

sts <- ch_data$sts
sts <- cbind(sts, extract(grid, sts[, c("lon", "lat")])[, -c(1:4)])

prec <- ch_data$data[dOe]
prec <- coredata(prec)

library(doParallel)

gridPcp(prec = prec,
        grid = grid,
        dyncovars = dyncovars,
        sts = sts,
        model_fun = learner_glm,
        dates = dOe,
        ncpu = 100,
        thres = NA,
        neibs = 20,
        coords = c('lon','lat'),
        crs = 'EPSG:4326',
        coords_as_preds = TRUE,
        dir_name = "ch_glm")

gridPcp(prec = prec,
        grid = grid,
        dyncovars = dyncovars,
        sts = sts,
        model_fun = learner_rf,
        dates = dOe,
        ncpu = 100,
        thres = NA,
        neibs = 20,
        coords = c('lon','lat'),
        crs = 'EPSG:4326',
        coords_as_preds = TRUE,
        dir_name = "ch_rf")

gridPcp(prec = prec,
        grid = grid,
        dyncovars = dyncovars,
        sts = sts,
        model_fun = learner_xgboost,
        dates = dOe,
        ncpu = 100,
        thres = NA,
        neibs = 20,
        coords = c('lon','lat'),
        crs = 'EPSG:4326',
        coords_as_preds = TRUE,
        dir_name = "ch_xgboost")

# to plot
# renv::install("ropensci/rnaturalearthhires")
shp_ch_data <- ne_countries(
  scale = 10,
  type = "countries",
  country = "Switzerland",
  returnclass = c("sf")
)
shp_ch_data <- vect(as(shp_ch_data, "Spatial"))

# to compare new grids
# rd <- rast("/mnt/climstor/meteoswiss/RhydchprobD/LV95/2014/RhydchprobD_grid.ens_ch01h.swiss.lv95_201407010000_201407310000.nc")
# rd_24 <- median(rd[[paste("RhydchprobD_realization=", 1:50,"_24", sep = "")]])
# rd_25 <- median(rd[[paste("RhydchprobD_realization=", 1:50,"_25", sep = "")]])
# rd_26 <- median(rd[[paste("RhydchprobD_realization=", 1:50,"_26", sep = "")]])
# rd_27 <- median(rd[[paste("RhydchprobD_realization=", 1:50,"_27", sep = "")]])
# rd_28 <- median(rd[[paste("RhydchprobD_realization=", 1:50,"_28", sep = "")]])
# rd_ev <- c(rd_24, rd_25, rd_26, rd_27, rd_28)
# rd_ev <- project(rd_ev, grid)
# 
# terra::writeCDF(rd_ev,
#                 "data/grid_prc_data/RhydchprobD_median_event.nc")

# precip fields
xe_glm <- terra::rast(dir("./pred_ch_glm/", full.names = TRUE))
xe_rf <- terra::rast(dir("./pred_ch_rf/", full.names = TRUE))
xe_xgboost <- terra::rast(dir("./pred_ch_xgboost/", full.names = TRUE))
rd_ev <- terra::rast("data/grid_prc_data/RhydchprobD_median_event.nc")
names(rd_ev) <- names(xe_glm)

# ### maps
xe_glm_df <- reshape2::melt(as.data.frame(mask(xe_glm, shp_ch_data), xy = TRUE), id.vars = c("x", "y"), value.name = "prec", variable.name = "date_time")
xe_glm_df$model <- "glm"
xe_rf_df <- reshape2::melt(as.data.frame(mask(xe_rf, shp_ch_data), xy = TRUE), id.vars = c("x", "y"), value.name = "prec", variable.name = "date_time")
xe_rf_df$model <- "rf"
xe_xgboost_df <- reshape2::melt(as.data.frame(mask(xe_xgboost, shp_ch_data), xy = TRUE), id.vars = c("x", "y"), value.name = "prec", variable.name = "date_time")
xe_xgboost_df$model <- "xgboost"
xe_rd_df <- reshape2::melt(as.data.frame(mask(rd_ev, shp_ch_data), xy = TRUE), id.vars = c("x", "y"), value.name = "prec", variable.name = "date_time")
xe_rd_df$model <- "RhydchprobD"

xe_all <- rbind(xe_glm_df, xe_rf_df, xe_xgboost_df, xe_rd_df)
xe_all$model <- factor(xe_all$model, levels = c("glm", "rf", "xgboost", "RhydchprobD"))
xe_all$date_time <- factor(xe_all$date_time)
xe_all$prec_m <- cut(xe_all$prec, breaks = c(0, 0.3, 2, 5, 10, 20, 35, 50, 70, 90, 120, 150), include.lowest = TRUE,
                     labels = c("[0,0.3]", "(0.3,2]" ,"(2,5]","(5,10]", "(10,20]", "(20,35]", "(35,50]", "(50,70]", "(70,90]", "(90,120]", "(120,150]"))
xe_all$prec_m <- factor(xe_all$prec_m)


ggplot() +
  geom_tile(data = xe_all, aes(x = x, y = y, fill = prec_m)) +
  geom_sf(data = sf::st_geometry(sf::st_as_sf(shp_ch_data)), fill = NA, color = "black") +
  scale_fill_manual(values = c("#ffffff", "#e8fbc2", "#c9ffcd", "#94f1b1", "#53bd9e", "#37a696", "#4495b4", "#356fb0", "#2a4f8c", "#331796", "#310146"),
                    drop = FALSE) +
  facet_grid(model~date_time, switch = "y") +
  theme_bw() +
  theme(legend.position ="bottom",
        axis.title =element_blank(),
        axis.text = element_blank(),
        axis.ticks =element_blank(),
        legend.key.height = unit(2, "mm"),
        plot.margin = grid::unit(c(-5, 0.5, -5, 0), "mm"),
        legend.key = element_rect(color = "black"),
        legend.margin = margin(c(-5, 0, .5, 0))) +
  guides(fill = guide_legend(label.position = "bottom", nrow = 1, title = "mm/day"))

ggsave("output/Fig_02-model-grid-comparison-RhydchprobD.pdf", width = 700, heigh = 500, scale = 3, units = "px")


### scatterplot
obs_prec_dOe <- data.frame(t(ch_data$data[dOe]))
colnames(obs_prec_dOe) <- dOe
obs_prec_dOe$ID <- rownames(obs_prec_dOe)
obs_prec_dOe <- reshape2::melt(obs_prec_dOe, id.vars = c("ID"), value.name = "obs", variable.name = "date_time")

rd_prec_dOe <- extract(rd_ev, ch_data$sts[, c("lon", "lat")])
rd_prec_dOe$ID <- ch_data$sts$ID
rd_prec_dOe <- reshape2::melt(rd_prec_dOe, id.vars = c("ID"), value.name = "model", variable.name = "date_time")
rd_prec_dOe$grid <- "RhydchprobD"

glm_prec_dOe <- extract(xe_glm, ch_data$sts[, c("lon", "lat")])
glm_prec_dOe$ID <- ch_data$sts$ID
glm_prec_dOe <- reshape2::melt(glm_prec_dOe, id.vars = c("ID"), value.name = "model", variable.name = "date_time")
glm_prec_dOe$grid <- "glm"

rf_prec_dOe <- extract(xe_rf, ch_data$sts[, c("lon", "lat")])
rf_prec_dOe$ID <- ch_data$sts$ID
rf_prec_dOe <- reshape2::melt(rf_prec_dOe, id.vars = c("ID"), value.name = "model", variable.name = "date_time")
rf_prec_dOe$grid <- "rf"

xgboost_prec_dOe <- extract(xe_xgboost, ch_data$sts[, c("lon", "lat")])
xgboost_prec_dOe$ID <- ch_data$sts$ID
xgboost_prec_dOe <- reshape2::melt(xgboost_prec_dOe, id.vars = c("ID"), value.name = "model", variable.name = "date_time")
xgboost_prec_dOe$grid <- "xgboost"

prec_dOe_all <- rbind(
  merge(obs_prec_dOe, rd_prec_dOe, by = c("ID", "date_time")),
  merge(obs_prec_dOe, glm_prec_dOe, by = c("ID", "date_time")),
  merge(obs_prec_dOe, rf_prec_dOe, by = c("ID", "date_time")),
  merge(obs_prec_dOe, xgboost_prec_dOe, by = c("ID", "date_time"))
)

prec_dOe_all$grid <- factor(prec_dOe_all$grid, levels = c("glm", "rf", "xgboost", "RhydchprobD"))

dr_bcc_dOe_all <- lapply(unique(prec_dOe_all$date_time),
                         function(x){
                           
                           df_t <- prec_dOe_all[prec_dOe_all$date_time == x, ]
                           
                           out <- lapply(unique(prec_dOe_all$grid),
                                         function(z){
                                           
                                           df_t_z <- df_t[df_t$grid == z, c("obs", "model")]
                                           colnames(df_t_z) <- c("obs", "mod")
                                           df_t_z <- df_t_z[complete.cases(df_t_z), ]
                                           data.frame(format(round(get_dr_bcc(df_t_z), digits = 2), nsmall = 2), grid = z)
                                           
                                         })
                           
                           data.frame(do.call(rbind, out), date_time = x)
                           
                         })
dr_bcc_dOe_all <- do.call(rbind, dr_bcc_dOe_all)

ggplot() +
  geom_point(data = prec_dOe_all, aes(x = obs, y = model), shape = 21) +
  geom_text(data = dr_bcc_dOe_all,
            aes(x = 15, y = 65, label = paste("dr  : ", dr)),
            size = 3, color = "gray30") +
  geom_text(data = dr_bcc_dOe_all,
            aes(x = 15, y = 57, label = paste("bcc: ", bcc)),
            size = 3, color = "gray30") +
  facet_grid(grid~date_time, switch = "y") +
  theme_bw()

ggsave("output/Fig_03-model-grid-comparison-scatterplot-RhydchprobD.pdf", width = 700, heigh = 500, scale = 3, units = "px")


# supplementary
sf_01 <- dyncovars[[1]]
names(sf_01) <- paste(names(sf_01), dOe, sep = "_")
sf_02 <- dyncovars[[2]]
names(sf_02) <- paste(names(sf_02), dOe, sep = "_")

pdf("output/FigS_01-ch.pdf", width = 9, height = 6.5)
plot(c(grid, sf_01, sf_02), axes = FALSE, legend = FALSE, mar = c(0, 0, 0, 0), nc = 5, maxnl = 18)
dev.off()
