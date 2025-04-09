rm(list = ls())

library(xts)
library(ggplot2)
library(reddPrec)
source("R/metrics.R")

# swiss example data
ch_data <- readRDS("data/ch_gf_data.RDS")

# swiss predictor data
ch_data_covs <- prcomp(scale(ch_data$covs[, -c(9:11)]))
ch_data_covs <- cbind(ch_data$sts, ch_data_covs$x[, 1:4])

prec <- zoo::coredata(ch_data$data)
sts <- ch_data_covs

# gap-filling methods
filled_glm <- gapFilling(
  prec, sts, model_fun = learner_glm,
  dates = time(ch_data$data),
  stmethod = "ratio", thres = NA, coords = c('lon','lat'),
  coords_as_preds = TRUE, crs = 'EPSG:4326', neibs = 15,
  window = 11, ncpu = 68)

filled_rf <- gapFilling(
  prec, sts, model_fun = learner_rf,
  dates = seq.Date(as.Date('2010-01-01'), as.Date('2015-12-31'), by = 'day'),
  stmethod = "ratio", thres = NA, coords = c('lon','lat'),
  coords_as_preds = TRUE, crs = 'EPSG:4326', neibs = 15,
  window = 11, ncpu = 68)

filled_xgboost <- gapFilling(
  prec, sts, model_fun = learner_xgboost,
  dates = seq.Date(as.Date('2010-01-01'), as.Date('2015-12-31'), by = 'day'),
  stmethod = "ratio", thres = NA, coords = c('lon','lat'),
  coords_as_preds = TRUE, crs = 'EPSG:4326', neibs = 15,
  window = 11, ncpu = 68)

# computing metrics
glm_stat <-
  lapply(unique(filled_glm$ID),
         function(x) {
           
           df_obs_mod <- filled_glm[filled_glm$ID == x, ]
           df_obs_mod <- df_obs_mod[, c("obs", "mod_pred")]
           df_obs_mod <- df_obs_mod[complete.cases(df_obs_mod), ]
           colnames(df_obs_mod) <- c("obs", "mod")
           
           out <- get_dr_bcc(df_obs_mod)
           out$ID <- x
           out
           
         })
glm_stat <- do.call(rbind, glm_stat)

rf_stat <-
  lapply(unique(filled_rf$ID),
         function(x) {
           
           df_obs_mod <- filled_rf[filled_glm$ID == x, ]
           df_obs_mod <- df_obs_mod[, c("obs", "mod_pred")]
           df_obs_mod <- df_obs_mod[complete.cases(df_obs_mod), ]
           colnames(df_obs_mod) <- c("obs", "mod")
           
           out <- get_dr_bcc(df_obs_mod)
           out$ID <- x
           out
           
         })
rf_stat <- do.call(rbind, rf_stat)

xgboost_stat <-
  lapply(unique(filled_xgboost$ID),
         function(x) {
           
           df_obs_mod <- filled_xgboost[filled_glm$ID == x, ]
           df_obs_mod <- df_obs_mod[, c("obs", "mod_pred")]
           df_obs_mod <- df_obs_mod[complete.cases(df_obs_mod), ]
           colnames(df_obs_mod) <- c("obs", "mod")
           
           out <- get_dr_bcc(df_obs_mod)
           out$ID <- x
           out
           
         })
xgboost_stat <- do.call(rbind, xgboost_stat)

# plot
df_stat <- rbind(data.frame(glm_stat, model = "glm"),
                 data.frame(rf_stat, model = "rf"),
                 data.frame(xgboost_stat, model = "xgboost"))

df_stat <- reshape2::melt(
  df_stat,
  id.vars = c("model", "ID"),
  variable.name = "metrics"
)

ggplot(df_stat,
       aes(x = model, y = value, colour = model)) +
  geom_point(position = position_jitterdodge()) +
  geom_boxplot() +
  scale_x_discrete(drop = FALSE) +
  theme_bw() +
  facet_wrap( ~factor(metrics,
                      c("bcc", "dr"))) +
  theme(legend.position = "none") + 
  scale_y_continuous(limits = c(0.5, 1))

ggsave("output/Fig_01-model-gapfilling-comparison.pdf", width = 666, heigh = 300, scale = 2.2, units = "px")

