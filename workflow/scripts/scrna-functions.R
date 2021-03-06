#!/usr/bin/env Rscript



function_pca_dimensions=function(Spatial_Data){
  
  pct <- Stdev(object = Spatial_Data, reduction = "pca") / sum(Stdev(object = Spatial_Data, reduction = "pca")) * 100
# Calculate cumulative percents for each PC
cum <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cum > 90 & pct < 5)[1]
co1
co2 <- sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > 0.1),  decreasing = T)[1] + 1 
# last point where change of % of variation is more than 0.1%.
co2
# Minimum of the two calculation
dimensionReduction <- min(co1, co2) 
  
}