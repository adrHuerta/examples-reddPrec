rm(list = ls())

library(terra)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyterra)
library(ggplot2)
library(reddPrec)
source("R/metrics.R")

# Valencia event daily data 2024-10-29 https://digital.csic.es/handle/10261/374208
daily_val <- read.csv("data/valencia_daily.csv")
daily_nc <- rast("data/grid_prc_data/valencia_daily.nc")

# topo data
topo_vars <- rast("data/topo_data/topo_vars_valencia.nc")
names(topo_vars) <- c('aspectcosine', 'aspectsine', 'dist2coast', 'dx', 'dxx', 'dy', 'dyy', 'eastness',
                      'elevation', 'latitude', 'longitude', 'northness', 'pcurv', 'roughness', 'slope', 'tcurv', 'tpi',
                      'tri', 'vrm')
topo_vars <- topo_vars[[c(11, 10, 9)]]
names(topo_vars)[1:2] <- c("lon", "lat")


# surface reflectance data
sr <- rast("data/modis_data/MODIS_Aqua_sur_refl_b01_b02_scaled.tif")
sr_r <- resample(sr, topo_vars[[1]])
sr_r <- focal(sr_r, 3, "modal", na.policy = "only")

# grid data
grid <- c(topo_vars, sr_r)
grid <- crop(grid, ext(-3.5, 0.75, 36.5, 41.5))

# sts data
sts_all <- daily_val[daily_val$lon > -4 & daily_val$lon < 1 & daily_val$lat > 36 & daily_val$lat < 42,]
sts_all <- sts_all[complete.cases(sts_all), c("id", "lon", "lat", "pre")]
colnames(sts_all)[1] <- c("ID")
sts <- cbind(sts_all, extract(grid[[-c(1:2)]], sts_all[, c("lon", "lat")]))
sts <- sts[, c("ID", "lon", "lat", "elevation", "sur_refl_b01", "sur_refl_b02")]

# pre data
prec <- matrix(sts_all[,"pre"], nrow = 1)
colnames(prec) <- sts_all$ID

library(doParallel)

gridPcp(prec = prec, grid = grid, sts = sts,
        model_fun = learner_glm,
        dates = as.Date('2024-10-29'),
        thres = NA, coords = c('lon','lat'), coords_as_preds = TRUE, 
        crs = 'EPSG:4326', neibs = 15, ncpu = 100, dir_name = "va_glm")

gridPcp(prec = prec, grid = grid, sts = sts,
        model_fun = learner_rf,
        dates = as.Date('2024-10-29'),
        thres = NA, coords = c('lon','lat'), coords_as_preds = TRUE, 
        crs = 'EPSG:4326', neibs = 15, ncpu = 100, dir_name = "va_rf")

gridPcp(prec = prec, grid = grid, sts = sts,
        model_fun = learner_xgboost,
        dates = as.Date('2024-10-29'),
        thres = NA, coords = c('lon','lat'), coords_as_preds = TRUE, 
        crs = 'EPSG:4326', neibs = 15, ncpu = 100, dir_name = "va_xgboost")

# shp
spain_provinces <- ne_states(country = "Spain", returnclass = "sf")

# grids
r <- terra::rast(c('./pred_va_glm/20241029.tif', './pred_va_rf/20241029.tif', './pred_va_xgboost/20241029.tif'))
names(r) <- c("glm", "rf", "xgboost")

# to compare
daily_nc_r <- project(daily_nc, r)
names(daily_nc_r) <- "CSIC"

# to plot
df_rs <- c(daily_nc_r, r)
df_rs <- mask(df_rs, spain_provinces)
df_rs <- as.data.frame(df_rs, xy = TRUE)
df_rs <- reshape2::melt(df_rs, id.vars = c("x", "y"), value.name = "prec", variable.name = "grid")
df_rs <- df_rs[complete.cases(df_rs), ]
df_rs$prec_cut <-  cut(df_rs$prec,
                       breaks = c(0, 1, 5, 50, 100, 200, 300, 400, 500, 600, 700, 802), include.lowest = TRUE)
df_rs$grid <- factor(df_rs$grid, levels = c("glm", "rf", "xgboost", "CSIC"))

spc <-
  ggplot() +
  geom_tile(data = df_rs, aes(x = x, y = y, fill = prec_cut)) +
  scale_fill_manual(values = c("#ffffff", "#e8fbc2", "#c9ffcd", "#94f1b1", "#53bd9e", "#37a696", "#4495b4", "#356fb0", "#2a4f8c", "#331796", "#310146"),
                    drop = FALSE) +
  facet_wrap(~grid, ncol = 4) +
  geom_spatvector(data = spain_provinces, fill = NA) +
  scale_y_continuous(limits = c(38.8, 40.25), expand = c(0, 0)) + 
  scale_x_continuous(limits = c(-1.75, 0), expand = c(0, 0)) +
  ylab("") + theme_bw() +
  theme(legend.position ="bottom",
        axis.title =element_blank(),
        axis.text = element_blank(),
        axis.ticks =element_blank(),
        legend.key.height = unit(2, "mm"),
        plot.margin = grid::unit(c(-5, 0.5, -5, 0), "mm"),
        legend.key = element_rect(color = "black"),
        legend.margin = margin(c(-5, 0, .5, 0))) +
  guides(fill = guide_legend(label.position = "bottom", nrow = 1, title = "mm/day")) +
  guides(fill = guide_legend(label.position = "bottom", nrow = 1, title = "mm/day"))


# scatterplot

df_all <- data.frame(extract(c(daily_nc_r, r), sts[, c("lon", "lat")])[, -1], obs = prec[1,])
df_all <- reshape2::melt(df_all, measure.vars = c("CSIC", "glm", "rf", "xgboost"), variable.name = "grid", value.name = "model")
df_all$grid <- factor(df_all$grid, levels = c("glm", "rf", "xgboost", "CSIC"))

df_all_dr_bcc <- lapply(unique(df_all$grid),
                        function(j){
                          
                          jdf <- df_all[df_all$grid == j, c("obs", "model")]
                          colnames(jdf) <- c("obs", "mod")
                          jdf <- jdf[complete.cases(jdf),]
                          
                          # as the CSIS does not deal with negative values i just convert all negatives to zero
                          if(any(jdf$mod < 0)){
                            
                            jdf <- jdf[-which(jdf$mod < 0), ]
                            
                          }
                          
                          data.frame(format(round(get_dr_bcc(jdf), digits = 2), nsmall = 2), grid = j)
                          
                        })
df_all_dr_bcc <- do.call(rbind, df_all_dr_bcc)

bpc <-
  ggplot() + 
  geom_point(data = df_all, aes(x = obs, y = model), shape = 21) +
  geom_text(data = df_all_dr_bcc,
            aes(x = 150, y = 800-10, label = paste("dr  : ", dr)),
            size = 3, color = "gray30") +
  geom_text(data = df_all_dr_bcc,
            aes(x = 150, y = 725-10, label = paste("bcc: ", bcc)),
            size = 3, color = "gray30") +
  facet_wrap(~ grid, ncol = 4) + theme_bw() + 
  scale_x_continuous(limits = c(0, 800)) + 
  scale_y_continuous(limits = c(0, 800))

ggpubr::ggarrange(spc, bpc, ncol = 1, nrow = 2, widths = c(1, 0.1), heights = c(1, .8),
                  align = 'v', labels = c("a)", "b)"),
                  font.label = list(size = 12, face = "plain"))
ggsave("output/Fig_03-model-grid-comparison-scatterplot-CSIC.pdf", width = 730, heigh = 500, scale = 3, units = "px")


# supplementary
pdf("output/FigS_02-va.pdf", width = 9, height = 3.5)
local({
  counter <- 0
  plot(grid, fun = function(x) {
    counter <<- counter + 1
    if (counter == 1) {
      points(sts[, c("lon", "lat")], col = "red", pch = 20)
    }
  }, axes = FALSE, legend = FALSE, mar = c(0, 0, 0, 0), nc = 5, maxnl = 18)
})
dev.off()

