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

function_gofunc=function(df,algorithm="weight01",statistics="ks",mapping="org.Hs.eg.db",ID="symbol",ontology = "BP") {
  

geneList <- df$p_val
names(geneList) <- df$gene
# Create topGOData object
GOdata <- new("topGOdata",
              ontology = ontology,
              allGenes = geneList,
              geneSelectionFun = function(x)x,
              annot = annFUN.org, mapping = mapping,ID=ID,nodeSize = 10)

resultsKS <- runTest(GOdata,algorithm = algorithm,statistic = statistics)


tab <- GenTable(GOdata, raw.p.value = resultsKS, topNodes = length(resultsKS@score), numChar = 120)

return(tab)
}


theme_cellsnake_classic <- function (base_size = 12, base_family = "Ubuntu") {
    theme_bw() %+replace% 
    theme(
      panel.grid.major  = element_line(color = "white"),
      #panel.background = element_rect(fill = "lightblue"),
      #panel.border = element_rect(color = "lightblue", fill = NA),
      #axis.line = element_line(color = "lightblue"),
      #axis.ticks = element_line(color = "lightblue"),
      #axis.text = element_text(color = "steelblue")
      )
}