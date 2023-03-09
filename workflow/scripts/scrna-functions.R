#!/usr/bin/env Rscript



function_pca_dimensions <- function(Spatial_Data) {
  pct <- Stdev(object = Spatial_Data, reduction = "pca") / sum(Stdev(object = Spatial_Data, reduction = "pca")) * 100
  # Calculate cumulative percents for each PC
  cum <- cumsum(pct)
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cum > 90 & pct < 5)[1]
  co1
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  # last point where change of % of variation is more than 0.1%.
  co2
  # Minimum of the two calculation
  dimensionReduction <- min(co1, co2)
}


c(
  "#c15035",
  "#79ce5d",
  "#8449c1",
  "#cbbc50",
  "#c75291",
  "#83c8b7",
  "#4d2f48",
  "#c59581",
  "#8a90c4",
  "#506138"
) -> palette10

palette30 <- c(
  "#d2a33f",
  "#6844cb",
  "#81dc3c",
  "#c54cdb",
  "#6cda74",
  "#de43ac",
  "#d2d840",
  "#422679",
  "#c4da81",
  "#974096",
  "#659a3d",
  "#6677cd",
  "#db622a",
  "#66d6a9",
  "#d83e46",
  "#8dd4d6",
  "#ce4b7a",
  "#3e602b",
  "#ce90d4",
  "#8b7137",
  "#373053",
  "#d0c9a1",
  "#7d3d62",
  "#608e79",
  "#88352a",
  "#779ec5",
  "#d3866a",
  "#324a45",
  "#c898a6",
  "#503029"
)


palette50 <- c(
  "#42566f",
  "#66e741",
  "#7744d9",
  "#a0e340",
  "#d24cdd",
  "#5bb734",
  "#3d2897",
  "#e1e73b",
  "#9638ae",
  "#5fe277",
  "#e344a0",
  "#579836",
  "#6a7de4",
  "#a5bb36",
  "#5d50a4",
  "#c2e778",
  "#632267",
  "#84df99",
  "#a64194",
  "#e0bd42",
  "#31245c",
  "#d6d58a",
  "#c783de",
  "#4c9e66",
  "#e2462a",
  "#62dabf",
  "#db445e",
  "#77cfe3",
  "#d7862e",
  "#5c8ec9",
  "#ac482c",
  "#649aa4",
  "#9e345e",
  "#c4dfc5",
  "#251e33",
  "#948d34",
  "#e08bbb",
  "#416126",
  "#bcb4da",
  "#2f3922",
  "#de986c",
  "#4f2330",
  "#d8b6ac",
  "#6c2d23",
  "#417061",
  "#c5787b",
  "#909a6f",
  "#876690",
  "#825b26",
  "#84675c"
)

function_gofunc <- function(df, algorithm = "weight01", statistics = "ks", mapping = "org.Hs.eg.db", ID = "symbol", ontology = "BP") {
  geneList <- df$p_val
  names(geneList) <- df$gene
  # Create topGOData object
  GOdata <- new("topGOdata",
    ontology = ontology,
    allGenes = geneList,
    geneSelectionFun = function(x) x,
    annot = annFUN.org, mapping = mapping, ID = ID, nodeSize = 10
  )

  resultsKS <- runTest(GOdata, algorithm = algorithm, statistic = statistics)


  tab <- GenTable(GOdata, raw.p.value = resultsKS, topNodes = length(resultsKS@score), numChar = 120)

  return(tab)
}


theme_cellsnake_classic <- function(base_size = 12, base_family = "Ubuntu") {
  theme_bw() %+replace%
    theme(
      panel.grid.major = element_line(color = "white"),
      # panel.background = element_rect(fill = "lightblue"),
      # panel.border = element_rect(color = "lightblue", fill = NA),
      # axis.line = element_line(color = "lightblue"),
      # axis.ticks = element_line(color = "lightblue"),
      # axis.text = element_text(color = "steelblue")
    )
}


function_color_palette <- function(n) {
  if (n > 10 && n <= 30) {
    return(palette30)
  }
  if (n > 30) {
    return(palette50)
  }

  if (n <= 10) {
    return(palette10)
  }


  randomcoloR::distinctColorPalette(n)
}