rm(list = ls())

library(xts)
library(zoo)
library(ggplot2)
library(reddPrec)
source("R/metrics2.R")  # it uses dr and mcc

# ----- Step 1: Define Data -----
set.seed(123)

# Swiss best precipitation data
ch_data <- readRDS(file.path(file.path(dirname(dirname(getwd())), "datasets", "observed_precipitation", "others", "datos_CH_1863-2020_v2.RDS")))
dates <- seq(as.Date("1960-01-01"), as.Date("2015-12-31"), by = "day")

# Issue with SOG (Soglio) station by SB
ch_data$data <- ch_data$data[, -match("SOG", colnames(ch_data$data))]
ch_data$xyz <- ch_data$xyz[-match("SOG", ch_data$xyz$ID), ]
all(ch_data$xyz$ID == colnames(ch_data$data))

precip_xts <- xts::xts(ch_data$data, dates)
precip_xts <- precip_xts["1960/2015"]

# Amount of gaps
100 * (sum(is.na(precip_xts))/sum(!is.na(precip_xts)))

# Filling by daily mean values
filling_by_daily_means <- function(time_serie) {
  
  daily_yearly_cycle <- unique(format(time(time_serie), "%m-%d"))
  daily_yearly_cycle_precip <- sapply(daily_yearly_cycle, 
                                      function(x) mean(time_serie[format(time(time_serie), "%m-%d") %in% x], na.rm = TRUE))
  daily_yearly_cycle_precip <- round(daily_yearly_cycle_precip, 1)
  
  dates_NAs <- time(time_serie[is.na(time_serie)])
  daily_yearly_cycle_precip_to_fill <- daily_yearly_cycle_precip[format(dates_NAs, "%m-%d")]
  
  zoo::coredata(time_serie[dates_NAs]) <- daily_yearly_cycle_precip_to_fill
  
  return(time_serie)
} 

precip_xts <- parallel::mclapply(precip_xts,
                                 filling_by_daily_means,
                                 mc.cores = 68)
# Gap-filled time series
precip_xts <- do.call(cbind, precip_xts)
100 * (sum(is.na(precip_xts))/sum(!is.na(precip_xts)))

# ----- Step 2: Define Synthetic Break Injection Parameters -----
affected_percent <- 0.50                                                 # 50% of stations will be affected
available_dates <- dates[as.numeric(format(dates, "%Y")) %in% 1970:2000] # period in which a break will be added
apply_to_past <- TRUE                                                    # TRUE to apply changes to data before the break date

# Injection parameters for different break types
bias_value <- 2          # Constant bias (mm)
variance_factor <- 1.5   # Multiplicative factor for variability
wet_increase_prob <- 0.3 # Probability for increasing wet-day frequency
dry_increase_prob <- 0.3 # Probability for decreasing wet-day frequency

# Create a copy of the original series to inject synthetic breaks
precip_xts_wb <- precip_xts

# Data frame to store break metadata
break_info <- data.frame(Station = character(),
                         BreakDate = as.Date(character()),
                         InjectionType = character(),
                         Parameter = numeric(),
                         stringsAsFactors = FALSE)

# ----- Step 3: Select Stations and Inject Synthetic Breaks with Multiple Types -----
# Randomly select the stations to be affected based on the specified percentage
all_stations <- colnames(precip_xts)
num_affected <- ceiling(length(all_stations) * affected_percent)
affected_stations <- sample(all_stations, num_affected)

for (station in affected_stations) {
  # Choose a random break date for this station
  break_date <- sample(available_dates, 1)
  # Randomly choose an injection type: bias, variance, or frequency
  injection_type <- sample(c("bias", "variance", "frequency"), 1)
  
  # Define index based on whether you want to apply changes to the past or future
  if (apply_to_past) {
    idx <- index(precip_xts_wb) <= break_date
  } else {
    idx <- index(precip_xts_wb) > break_date
  }
  
  # Apply the selected injection type
  if (injection_type == "bias") {
    precip_xts_wb[idx, station] <- precip_xts_wb[idx, station] + bias_value
    param_used <- bias_value
    length_years <- length(unique(format(time(precip_xts_wb[idx, station]), "%Y")))
    
  } else if (injection_type == "variance") {
    precip_xts_wb[idx, station] <- precip_xts_wb[idx, station] * variance_factor
    param_used <- variance_factor
    length_years <- length(unique(format(time(precip_xts_wb[idx, station]), "%Y")))
    
  } else if (injection_type == "frequency") {
    # For frequency, decide whether to increase or decrease wet-day frequency
    if (runif(1) < 1) {# force to always dry_to_wet (0.5)
      # Increase wet-day frequency: convert some dry days (zeros) to wet days
      dry_idx <- which(precip_xts_wb[idx, station] == 0)
      new_vals <- rgamma(length(dry_idx), shape = 2, scale = 0.5)
      flip <- runif(length(dry_idx)) < wet_increase_prob
      if (any(flip)) {
        precip_xts_wb[idx, station][dry_idx[flip]] <- new_vals[flip]
      }
      param_used <- wet_increase_prob
      length_years <- length(unique(format(time(precip_xts_wb[idx, station]), "%Y")))
      injection_type <- paste(injection_type, "_", "dry_wet", sep = "")
      
    } else {
      # Decrease wet-day frequency: convert some wet days to dry days (set to zero)
      wet_idx <- which(precip_xts_wb[idx, station] > 0)
      flip <- runif(length(wet_idx)) < dry_increase_prob
      if (any(flip)) {
        precip_xts_wb[idx, station][wet_idx[flip]] <- 0
      }
      param_used <- dry_increase_prob
      length_years <- length(unique(format(time(precip_xts_wb[idx, station]), "%Y")))
      injection_type <- paste(injection_type, "_", "wet_dry", sep = "")
      
      
    }
  }
  
  
  # Record break metadata for the station
  break_info <- rbind(break_info,
                      data.frame(Station = station,
                                 BreakDate = break_date,
                                 BreakYears = length_years,
                                 InjectionType = injection_type,
                                 Parameter = param_used,
                                 stringsAsFactors = FALSE))
}

break_info

# # ----- Step 4: Visualize the Synthetic Breaks for Affected Stations -----
# for (station in unique(break_info$Station)) {
#   df_orig <- data.frame(Date = index(precip_xts), Precip = coredata(precip_xts)[, station])
#   df_synth <- data.frame(Date = index(precip_xts_wb), Precip = coredata(precip_xts_wb)[, station])
#   
#   p <- ggplot() +
#     geom_line(data = df_orig, aes(x = Date, y = Precip), color = "blue") +
#     geom_line(data = df_synth, aes(x = Date, y = Precip), color = "red", alpha = 0.7) +
#     labs(title = paste("Station:", station,
#                        "\nOriginal (Blue) vs Synthetic Break (Red)",
#                        "\nBreak Date:", break_info$BreakDate[break_info$Station == station],
#                        "\nBreak Years:", break_info$BreakYears[break_info$Station == station],
#                        "\nType:", break_info$InjectionType[break_info$Station == station]),
#          x = "Date", y = "Precipitation (mm)") +
#     theme_minimal()
#   
#   print(p)
# }

# ----- Step 5: Homogenization of corrupted Swiss data-----
precip_xts_wb_hmg <- hmg_Ts(prec = precip_xts_wb,
                            sts = ch_data$xyz,
                            neibs_min = 2,
                            neibs_max = 12,
                            cor_neibs = 0.5,
                            wet_day = -1,
                            perc_break = 22,
                            apply_qc = 0.5,
                            mm_apply_qc = 0.1,
                            window_c = 15,
                            ncpu = 68)

## comparison of results (only added breaks): detection
detection_res <- cbind(precip_xts_wb_hmg$det_data[match(break_info$Station, precip_xts_wb_hmg$det_data$ID), ], break_info)
detection_res$BreakYear <- as.numeric(format(break_info$BreakDate, "%Y"))
detection_res <- detection_res[complete.cases(detection_res), ]
BDA_res <- calculate_BDA(true_breaks = detection_res$BreakYear, detected_breaks = detection_res$year_of_break)
TA_res <- calculate_TA(true_breaks = detection_res$BreakYear, detected_breaks = detection_res$year_of_break)
detection_res <- data.frame(BDA = BDA_res, TA = TA_res, ID = NA)

## comparison of results (only added breaks): adjustment
adjusment_res <- lapply(break_info$Station,
                        function(idx){
                          
                          station_df <- break_info[break_info$Station == idx, ]
                          station_df_year <- format(station_df$BreakDate, "%Y")
                          
                          adj_synthetic_break_ts <- precip_xts_wb_hmg$adj_data[, idx][paste("/", station_df_year, sep = "")]
                          org_data_ts <- precip_xts[, idx][paste("/", station_df_year, sep = "")]
                          
                          res <- get_dr_bcc(data.frame(obs = as.numeric(org_data_ts), mod = as.numeric(adj_synthetic_break_ts)))
                          res$ID <- idx
                          
                          res
                          
                        })
adjusment_res <- do.call(rbind, adjusment_res)

## comparison of results (only added breaks): trends
trend_PRCPTOT_res <- lapply(break_info$Station,
                            function(idx){
                              
                              station_df <- break_info[break_info$Station == idx, ]
                              station_df_year <- format(station_df$BreakDate, "%Y")
                              
                              adj_synthetic_break_ts <- precip_xts_wb_hmg$adj_data[, idx][paste("/", station_df_year, sep = "")]
                              adj_synthetic_break_ts_prcptot <- xts::apply.yearly(adj_synthetic_break_ts, sum, na.rm = TRUE)
                              
                              org_data_ts <- precip_xts[, idx][paste("/", station_df_year, sep = "")]
                              org_data_ts_prcptot <- xts::apply.yearly(org_data_ts, sum, na.rm = TRUE)
                              
                              trend_obs <- lm(zoo::coredata(org_data_ts_prcptot) ~ as.numeric(format(time(org_data_ts_prcptot), "%Y")))
                              trend_adj_synthetic <- lm(zoo::coredata(adj_synthetic_break_ts_prcptot) ~ as.numeric(format(time(adj_synthetic_break_ts_prcptot), "%Y")))
                              
                              data.frame(obs = as.numeric(coef(trend_obs)[2]),
                                         mod = as.numeric(coef(trend_adj_synthetic)[2]))
                              
                            })
trend_PRCPTOT_res <- do.call(rbind, trend_PRCPTOT_res)

trend_r1mm_res <- lapply(break_info$Station,
                         function(idx){
                           
                           station_df <- break_info[break_info$Station == idx, ]
                           station_df_year <- format(station_df$BreakDate, "%Y")
                           
                           adj_synthetic_break_ts <- precip_xts_wb_hmg$adj_data[, idx][paste("/", station_df_year, sep = "")]
                           adj_synthetic_break_ts_prcptot <- xts::apply.yearly(adj_synthetic_break_ts,  function(x) sum(x > 0.1, na.rm = TRUE))
                           
                           org_data_ts <- precip_xts[, idx][paste("/", station_df_year, sep = "")]
                           org_data_ts_prcptot <- xts::apply.yearly(org_data_ts,  function(x) sum(x > 0.1, na.rm = TRUE))
                           
                           trend_obs <- lm(zoo::coredata(org_data_ts_prcptot) ~ as.numeric(format(time(org_data_ts_prcptot), "%Y")))
                           trend_adj_synthetic <- lm(zoo::coredata(adj_synthetic_break_ts_prcptot) ~ as.numeric(format(time(adj_synthetic_break_ts_prcptot), "%Y")))
                           
                           data.frame(obs = as.numeric(coef(trend_obs)[2]),
                                      mod = as.numeric(coef(trend_adj_synthetic)[2]))
                           
                         })
trend_r1mm_res <- do.call(rbind, trend_r1mm_res)

trend_PRCPTOT_res <- get_dr_bcc(trend_PRCPTOT_res)
colnames(trend_PRCPTOT_res) <- c("dr-slope", "mcc-slope")
trend_PRCPTOT_res$ID <- NA
trend_r1mm_res <- get_dr_bcc(trend_r1mm_res)
colnames(trend_r1mm_res) <- c("dr-slope", "mcc-slope")
trend_r1mm_res$ID <- NA

## comparison of results (only added breaks): plot

to_plot <- rbind(data.frame(reshape2::melt(detection_res, measure.vars = c("BDA", "TA"), variable.name = "Metric"), type = "Detection"),
                 data.frame(reshape2::melt(adjusment_res, measure.vars = c("dr", "mcc"), variable.name = "Metric"), type = "Adjustment"),
                 data.frame(reshape2::melt(trend_PRCPTOT_res, measure.vars = c("dr-slope", "mcc-slope"), variable.name = "Metric"), type = "Trend-PRCPTOT"),
                 data.frame(reshape2::melt(trend_r1mm_res, measure.vars = c("dr-slope", "mcc-slope"), variable.name = "Metric"), type = "Trend-R1mm"))
to_plot$type <- factor(to_plot$type, levels = c("Detection", "Adjustment", "Trend-PRCPTOT", "Trend-R1mm"))

to_plot_mean <- aggregate(value ~ Metric + type, data = to_plot, FUN = function(x) round(mean(x), 2))

ggplot() +
  geom_boxplot(data = to_plot, aes(y = value, x = Metric), outliers = FALSE) + 
  facet_grid(~type, scales = "free_x") + 
  geom_text(data = to_plot_mean, aes(y = value + 0.04, x = Metric, label = value), size = 3)+
  scale_y_continuous(limits = c(0.45, 1.05)) +
  xlab("") + ylab("") +
  theme_bw()

ggsave("output/Fig_07-hmg-corruptedSwiss-and-Swiss-metrics-only-added-breaks.pdf", width = 600, heigh = 250, scale = 3, units = "px")


## comparison of results (whole time series): trends
trend_PRCPTOT_all <- lapply(ch_data$xyz$ID,
                            function(idx){
                              
                              adj_synthetic_break_ts <- precip_xts_wb_hmg$adj_data[, idx]
                              adj_synthetic_break_ts_prcptot <- xts::apply.yearly(adj_synthetic_break_ts, sum, na.rm = TRUE)
                              
                              org_data_ts <- precip_xts[, idx]
                              org_data_ts_prcptot <- xts::apply.yearly(org_data_ts, sum, na.rm = TRUE)
                              
                              trend_obs <- lm(zoo::coredata(org_data_ts_prcptot) ~ as.numeric(format(time(org_data_ts_prcptot), "%Y")))
                              trend_adj_synthetic <- lm(zoo::coredata(adj_synthetic_break_ts_prcptot) ~ as.numeric(format(time(adj_synthetic_break_ts_prcptot), "%Y")))
                              
                              data.frame(obs = as.numeric(coef(trend_obs)[2]),
                                         mod = as.numeric(coef(trend_adj_synthetic)[2]),
                                         ID = idx)
                              
                            })
trend_PRCPTOT_all <- do.call(rbind, trend_PRCPTOT_all)
trend_PRCPTOT_all$Corrupted <- "No"
trend_PRCPTOT_all[match(break_info$Station, trend_PRCPTOT_all$ID), ]$Corrupted <- "Yes"

# lattice::xyplot(org_data_ts[, idx], type = "p", cex = .1)
# lattice::xyplot(precip_xts_wb[, idx], type = "p", cex = .1)
# lattice::xyplot(precip_xts_wb_hmg$adj_data[, idx], type = "p", cex = .1)

trend_r1mm_all <- lapply(ch_data$xyz$ID,
                         function(idx){
                           
                           adj_synthetic_break_ts <- precip_xts_wb_hmg$adj_data[, idx]
                           adj_synthetic_break_ts_prcptot <- xts::apply.yearly(adj_synthetic_break_ts,  function(x) sum(x > 0.1, na.rm = TRUE))
                           
                           org_data_ts <- precip_xts[, idx]
                           org_data_ts_prcptot <- xts::apply.yearly(org_data_ts,  function(x) sum(x > 0.1, na.rm = TRUE))
                           
                           trend_obs <- lm(zoo::coredata(org_data_ts_prcptot) ~ as.numeric(format(time(org_data_ts_prcptot), "%Y")))
                           trend_adj_synthetic <- lm(zoo::coredata(adj_synthetic_break_ts_prcptot) ~ as.numeric(format(time(adj_synthetic_break_ts_prcptot), "%Y")))
                           
                           data.frame(obs = as.numeric(coef(trend_obs)[2]),
                                      mod = as.numeric(coef(trend_adj_synthetic)[2]),
                                      ID = idx)
                           
                         })
trend_r1mm_all <- do.call(rbind, trend_r1mm_all)
trend_r1mm_all$Corrupted <- "No"
trend_r1mm_all[match(break_info$Station, trend_r1mm_all$ID), ]$Corrupted <- "Yes"

## comparison of results (whole time series): trends - plot
trend_all <- rbind(data.frame(trend_PRCPTOT_all, index = "PRCPTOT"),
                   data.frame(trend_r1mm_all, index = "R1mm"))

trend_all_dr_bcc <- rbind(
  data.frame(get_dr_bcc(trend_all[trend_all$index == "PRCPTOT",]), index = "PRCPTOT", Corrupted = "ALL"),
  data.frame(get_dr_bcc(trend_all[trend_all$index == "PRCPTOT" & trend_all$Corrupted == "No",]), index = "PRCPTOT", Corrupted = "No"),
  data.frame(get_dr_bcc(trend_all[trend_all$index == "R1mm",]), index = "R1mm", Corrupted = "ALL"),
  data.frame(get_dr_bcc(trend_all[trend_all$index == "R1mm" & trend_all$Corrupted == "No",]), index = "R1mm", Corrupted = "No")
)

trend_all_dr_bcc$Corrupted <- factor(trend_all_dr_bcc$Corrupted, levels = c("No", "Yes", "ALL"))
trend_all_dr_bcc <- transform(trend_all_dr_bcc,
                              label = paste(" dr:", round(dr, 2), " /", " mcc:", round(mcc, 2), sep = ""))
trend_all_dr_bcc$mod <- c(-4, -5, -1.15, -1.35)
trend_all_dr_bcc$obs <- c(1.95, 1.95, 0.45, 0.47)

precip_xts_wb_hmg_trend_all_plt <-
  ggplot() +
  geom_point(data = trend_all, aes(x = obs, y = mod, color = Corrupted), size = 2) +
  scale_color_manual(values = c("#3C9AB2", "#F22300", "black"), na.value = "black", breaks = c("No", "Yes")) +
  guides(color = guide_legend(override.aes = list(shape = 16, size = 3))) + 
  xlab("Raw") + ylab("Homogenized") +
  facet_wrap(~ index, scales = "free") +
  theme_bw() + 
  theme(legend.position="bottom",
        legend.margin=margin(0,0,0,0),
        legend.box.spacing = unit(0, "pt")) +
  geom_text(data = trend_all_dr_bcc, aes(x = obs, y = mod, label = label, color = Corrupted))


# ----- Step 6: Homogenization of original Swiss data -----
precip_xts_hmg <- hmg_Ts(prec = precip_xts,
                         sts = ch_data$xyz,
                         neibs_min = 2,
                         neibs_max = 12,
                         perc_break = 22,
                         apply_qc = 0.5,
                         mm_apply_qc = 0.1,
                         window_c = 15,
                         ncpu = 68)

## comparison of results (whole time series): trends
trend_PRCPTOT_all <- lapply(ch_data$xyz$ID,
                            function(idx){
                              
                              adj_synthetic_break_ts <- precip_xts_hmg$adj_data[, idx]
                              adj_synthetic_break_ts_prcptot <- xts::apply.yearly(adj_synthetic_break_ts, sum, na.rm = TRUE)
                              
                              org_data_ts <- precip_xts[, idx]
                              org_data_ts_prcptot <- xts::apply.yearly(org_data_ts, sum, na.rm = TRUE)
                              
                              trend_obs <- lm(zoo::coredata(org_data_ts_prcptot) ~ as.numeric(format(time(org_data_ts_prcptot), "%Y")))
                              trend_adj_synthetic <- lm(zoo::coredata(adj_synthetic_break_ts_prcptot) ~ as.numeric(format(time(adj_synthetic_break_ts_prcptot), "%Y")))
                              
                              data.frame(obs = as.numeric(coef(trend_obs)[2]),
                                         mod = as.numeric(coef(trend_adj_synthetic)[2]),
                                         ID = idx)
                              
                            })
trend_PRCPTOT_all <- do.call(rbind, trend_PRCPTOT_all)

trend_r1mm_all <- lapply(ch_data$xyz$ID,
                         function(idx){
                           
                           adj_synthetic_break_ts <- precip_xts_hmg$adj_data[, idx]
                           adj_synthetic_break_ts_prcptot <- xts::apply.yearly(adj_synthetic_break_ts,  function(x) sum(x > 0.1, na.rm = TRUE))
                           
                           org_data_ts <- precip_xts[, idx]
                           org_data_ts_prcptot <- xts::apply.yearly(org_data_ts,  function(x) sum(x > 0.1, na.rm = TRUE))
                           
                           trend_obs <- lm(zoo::coredata(org_data_ts_prcptot) ~ as.numeric(format(time(org_data_ts_prcptot), "%Y")))
                           trend_adj_synthetic <- lm(zoo::coredata(adj_synthetic_break_ts_prcptot) ~ as.numeric(format(time(adj_synthetic_break_ts_prcptot), "%Y")))
                           
                           data.frame(obs = as.numeric(coef(trend_obs)[2]),
                                      mod = as.numeric(coef(trend_adj_synthetic)[2]),
                                      ID = idx)
                           
                         })
trend_r1mm_all <- do.call(rbind, trend_r1mm_all)

## comparison of results (whole time series): trends -plots
trend_all <- rbind(data.frame(trend_PRCPTOT_all, index = "PRCPTOT"),
                   data.frame(trend_r1mm_all, index = "R1mm"))

trend_all_dr_bcc <- rbind(
  data.frame(get_dr_bcc(trend_all[trend_all$index == "PRCPTOT",]), index = "PRCPTOT"),
  data.frame(get_dr_bcc(trend_all[trend_all$index == "R1mm",]), index = "R1mm")
)

trend_all_dr_bcc <- transform(trend_all_dr_bcc,
                              label = paste(" dr:", round(dr, 2), " /", " mcc:", round(mcc, 2), sep = ""))
trend_all_dr_bcc$mod <- c(-2, -.25)
trend_all_dr_bcc$obs <- c(2, 0.47)

precip_xts_trend_all_plt <-
  ggplot() +
  geom_point(data = trend_all, aes(x = obs, y = mod)) +
  facet_wrap(~ index, scales = "free") +
  theme_bw() +
  xlab("Raw") + ylab("Homogenized") +
  facet_wrap(~ index, scales = "free") +
  theme_bw()  +
  geom_text(data = trend_all_dr_bcc, aes(x = obs, y = mod, label = label))


ggpubr::ggarrange(precip_xts_wb_hmg_trend_all_plt,
                  precip_xts_trend_all_plt,
                  ncol = 1, nrow = 2, widths = c(1, 1), heights = c(1, .75),
                  align = 'v', labels = c("a)", "b)"),
                  font.label = list(size = 12, face = "plain"))
ggsave("output/Fig_08-hmg-corruptedSwiss-and-Swiss-trends.pdf", width = 650, heigh = 500, scale = 3, units = "px")
