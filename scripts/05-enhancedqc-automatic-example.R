rm(list = ls())

library(xts)
library(ggplot2)
library(rnaturalearth)
library(giscoR)
library(rnaturalearthdata)
library(reddPrec)

# ensuring 1 core is used (in my local I had to do it, check your own local)
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

# Switzerland (this data is private)
ch_data <- readRDS(file.path(file.path(dirname(dirname(getwd())), "datasets", "observed_precipitation", "others", "datos_CH_1863-2020_v2.RDS")))
ch_data$data <- xts::xts(ch_data$data, seq(as.Date("1960-01-01"), as.Date("2015-12-31"), by = "day"))

# automatic enhanced QC
library(doParallel)
ch_eqc <- eqc_Ts(prec = ch_data$data, ncpu = 50)

# plotting
all(ch_data$ID == rownames(ch_eqc))
ch_eqc <- cbind(ch_data$xyz[, c("LON", "LAT")], ch_eqc[, -1])
ch_eqc <- reshape2::melt(ch_eqc, c("LON", "LAT"), variable.name = "EnhancedQC", value.name = "Level")
ch_eqc$EnhancedQC <- factor(ch_eqc$EnhancedQC, levels = c("Truncation", "Small_Gaps", "Weekly_Cycle", "Precision_Rounding"))

dummy_eqc <- data.frame(LON = NA, LAT = NA, EnhancedQC = "Truncation", Level = 2)
ch_eqc <- rbind(ch_eqc, dummy_eqc)
ch_eqc$Level <- factor(ch_eqc$Level, levels = c(0, 1, 2), labels = c("0", "1", "2"))

# renv::install("ropensci/rnaturalearthhires")
shp_ch_data <- ne_countries(
  scale = 10,
  type = "countries",
  country = "Switzerland",
  returnclass = c("sf")
)

ch_plot <-
  ggplot() + 
  geom_sf(data = shp_ch_data, fill = NA, color = "black") +
  geom_point(data = ch_eqc, aes(x = LON, y = LAT, color = Level, size = Level)) + 
  scale_size_manual(values = c("0"=0.75,"1"=0.5, "2"=0.25)) +
  scale_color_manual(values = c("0"="skyblue","1"="lightgreen", "2"="red"), drop = FALSE) +
  facet_grid(~EnhancedQC) + 
  theme_bw() +
  theme(legend.position = "none",
        axis.title =element_blank(),
        axis.text = element_blank(),
        axis.ticks =element_blank())


# Spain (Aragon)
esp_data <- readRDS(file.path(file.path(dirname(dirname(getwd())), "datasets", "observed_precipitation", "others", "datos_ESParagon_1950_2020_v3.RDS")))
esp_data$data <- xts::xts(esp_data$data, seq(as.Date("1960-01-01"), as.Date("2015-12-31"), by = "day"))

# filtering stations based on Huerta et al. 2024
stations_longest_01 <- data.table::as.data.table(esp_data$data)
stations_longest_01$index <- format(stations_longest_01$index, "%m-%d")
stations_longest_01 <- stations_longest_01[, lapply(.SD, function(x) sum(!is.na(x))), by = index]
stations_longest_01 <- stations_longest_01[!(index == "02-29"), ]
stations_longest_01 <- stations_longest_01[, lapply(.SD, function(x) sum(x >= 10))]
stations_longest_01 <- unlist(stations_longest_01)[unlist(stations_longest_01) >= 365]

stations_longest_02 <- data.table::as.data.table(esp_data$data)
stations_longest_02$index <- format(stations_longest_02$index, "%Y")
stations_longest_02 <- stations_longest_02[, lapply(.SD, function(x) ifelse(sum(!is.na(x)) >= 300, 1, 0)), by = index]
stations_longest_02 <- stations_longest_02[, lapply(.SD, function(x){
  
  res <- rle(x)$lengths[which(rle(x)$values %in% 1)]
  any(res[res >= 5])
  
})]
stations_longest_02 <- unlist(stations_longest_02)[unlist(stations_longest_02) == TRUE]
stations_longest <- intersect(names(stations_longest_01), names(stations_longest_02))

esp_data$xyz$SIZE <- 0
esp_data$xyz[match(stations_longest, esp_data$xyz$ID), ]$SIZE <- 1

esp_data <- list(data = esp_data$data[, esp_data$xyz[esp_data$xyz$SIZE == 1, ]$ID],
                 xyz = esp_data$xyz[esp_data$xyz$SIZE == 1, ])
all(colnames(esp_data$data) == esp_data$xyz$ID)

# automatic enhanced QC
esp_eqc <- eqc_Ts(prec = esp_data$data, ncpu = 100)

# plotting
all(esp_data$xyz$ID == rownames(esp_eqc))
esp_eqc <- cbind(esp_data$xyz[, c("LON", "LAT")], esp_eqc[, -1])
esp_eqc <-reshape2::melt(esp_eqc, c("LON", "LAT"), variable.name = "EnhancedQC", value.name = "Level")
esp_eqc$EnhancedQC <- factor(esp_eqc$EnhancedQC, levels = c("Truncation", "Small_Gaps", "Weekly_Cycle", "Precision_Rounding"))

dummy_eqc <- data.frame(LON = NA, LAT = NA, EnhancedQC = "Truncation", Level = 2)
esp_eqc <- rbind(esp_eqc, dummy_eqc)
esp_eqc$Level <- factor(esp_eqc$Level, levels = c(0, 1, 2), labels = c("0", "1", "2"))

shp_esp_data <- gisco_get_nuts(country = "ES",
                               nuts_level = 2,
                               resolution = "3")

esp_plot <-
  ggplot() + 
  geom_sf(data = shp_esp_data, fill = NA, color = "black") +
  geom_point(data = esp_eqc, aes(x = LON, y = LAT, color = Level, size = Level)) + 
  scale_color_manual(values = c("0"="skyblue","1"="lightgreen", "2"="red"), drop = FALSE) +
  scale_size_manual(values = c("0"=.75,"1"=0.5, "2"=0.25)) +
  facet_grid(~EnhancedQC) + 
  theme_bw() +
  scale_x_continuous(limits = c(min(esp_eqc$LON, na.rm = TRUE), max(esp_eqc$LON, na.rm = TRUE))) +
  scale_y_continuous(limits = c(min(esp_eqc$LAT, na.rm = TRUE), max(esp_eqc$LAT, na.rm = TRUE))) +
  theme(legend.position = "bottom",
        axis.title =element_blank(),
        axis.text = element_blank(),
        axis.ticks =element_blank(),
        strip.text.x = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size = 2)))

ggpubr::ggarrange(ch_plot, esp_plot, ncol = 1, nrow = 2, widths = c(1, 1), heights = c(.91, 1.65),
                  align = 'v', labels = c(NA, NA),
                  font.label = list(size = 12, face = "plain"))

ggsave("output/Fig_06-enhancedqc-comparison-switzerland-spain(aragon).pdf", width = 630, heigh = 350, scale = 3, units = "px")
