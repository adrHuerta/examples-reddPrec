rm(list = ls())

library(xts)
library(reddPrec) 

ts_data <- readRDS("data/enhanced-qc-plots.RDS")

# a really good precipitation time series
pdf("output/Fig_05-enhancedqc-plot-verygood.pdf", width = 8, height = 7)
eqc_Plot(xts_obj = ts_data[, 1])
dev.off()
eqc_Truncation(xts_obj = ts_data[, 1])
eqc_SmallGaps(xts_obj = ts_data[, 1])
eqc_WeeklyCycle(xts_obj = ts_data[, 1])
eqc_PrecisionRounding(xts_obj = ts_data[, 1])

# not really good in this case
pdf("output/Fig_05-enhancedqc-plot-notverygood.pdf", width = 8, height = 7)
eqc_Plot(xts_obj = ts_data[, 2])
dev.off()
eqc_Truncation(xts_obj = ts_data[, 2])
eqc_SmallGaps(xts_obj = ts_data[, 2])
eqc_WeeklyCycle(xts_obj = ts_data[, 2])
eqc_PrecisionRounding(xts_obj = ts_data[, 2])