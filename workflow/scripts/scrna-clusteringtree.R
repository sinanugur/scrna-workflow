#!/usr/bin/env Rscript


option_list <- list(
  optparse::make_option(c("--scale.factor"),
    type = "integer", default = 10000,
    help = "Scale factor [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--nfeatures"),
    type = "integer", default = 2000,
    help = "Highly variable features [default= %default]", metavar = "integer"
  ),

        optparse::make_option(c("--variable.selection.method"),
    type = "character", default = "vst",
    help = "Find variable features selection method [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--rds"),
    type = "character", default = NULL,
    help = "RAW rds file of a Seurat object", metavar = "character"
  ),
  optparse::make_option(c("--normalization.method"),
    type = "character", default = "LogNormalize",
    help = "Normalization method[default= %default]", metavar = "character"
  ),

  optparse::make_option(c("--integration"), action = "store_true", default = FALSE),
  optparse::make_option(c("--clplot"),
    type = "character", default = "clustree.pdf",
    help = "Output clustree file name", metavar = "character"
  ),
  optparse::make_option(c("--jeplot"),
    type = "character", default = "jackandelbow.pdf",
    help = "Output jack and elbow file name", metavar = "character"
  ),
  optparse::make_option(c("--hvfplot"),
    type = "character", default = "variable-features.pdf",
    help = "Variable features file name", metavar = "character"
  ),
  optparse::make_option(c("--heplot"),
    type = "character", default = "dimheatmap.pdf",
    help = "Dim heatmap plot file name", metavar = "character"
  )
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

if (is.null(opt$rds)) {
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file and sampleid)", call. = FALSE)
}

require(tidyverse)
require(optparse)
require(Seurat)
require(clustree)
try({source("workflow/scripts/scrna-functions.R")},silent=TRUE)
try({source(paste0(system("python -c 'import os; import cellsnake; print(os.path.dirname(cellsnake.__file__))'", intern = TRUE),"/scrna/workflow/scripts/scrna-functions.R"))},silent=TRUE)

scrna <- readRDS(file = opt$rds)

if(isFALSE(opt$integration)) {


scrna <- NormalizeData(scrna, normalization.method = opt$normalization.method, scale.factor = opt$scale.factor)
scrna <- FindVariableFeatures(scrna, selection.method = opt$variable.selection.method, nfeatures = opt$nfeatures)
} else {

try({DefaultAssay(scrna) <- "integrated"}) #for now only for Seurat, Harmony will come

}


#all.genes <- rownames(scrna) memory requirements can be large if using all genes
not.all.genes <- VariableFeatures(scrna) #only variable features

scrna <- ScaleData(scrna, features = not.all.genes)
scrna <- RunPCA(scrna, features = not.all.genes)


if(isFALSE(opt$integration)) {
# output.dir=paste0("results/",opt$sampleid,"/technicals/")
# dir.create(output.dir,recursive = T)

# Identify the 10 most highly variable genes
top10 <- head(not.all.genes, 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(scrna)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)



# ggsave(paste0(output.dir,"highly-variable-features.pdf"), plot2 ,width = 8,height = 9)
ggsave(opt$hvfplot, plot2, width = 8, height = 9)

DimHeatmap(scrna, dims = 1:15, cells = 500, balanced = TRUE, fast = FALSE)
# ggsave(paste0(output.dir,"DimHeatMap_plot.pdf") ,width = 8,height = 15)
ggsave(opt$heplot, width = 8, height = 15)


scrna <- JackStraw(scrna, num.replicate = 20, dims = 50,verbose=FALSE)
scrna <- ScoreJackStraw(scrna, dims = 1:50)
plot1 <- JackStrawPlot(scrna, dims = 1:50)
plot2 <- ElbowPlot(scrna, ndims = 50)
# ggsave(paste0(output.dir,"JackandElbow_plot.pdf"), plot1 + plot2,width = 13,height = 5)
ggsave(opt$jeplot, plot1 + plot2, width = 13, height = 5)

}

if(isFALSE(opt$integration)) {

 resolution=seq(0.1, 2.5, 0.1)

} else {

resolution=seq(0.1, 1.5, 0.1)

  }
  


dimensionReduction <- function_pca_dimensions(scrna)
scrna <- FindNeighbors(scrna, dims = 1:dimensionReduction)
scrna <- FindClusters(scrna, resolution = resolution)

clustree(scrna) -> p1

# output.dir=paste0("results/",opt$sampleid,"/clusteringTree/")
# dir.create(output.dir,recursive = T)

# ggsave(paste0(output.dir,"/clusteringTree-",opt$sampleid,".pdf"),p1,width=8,height=15)
ggsave(opt$clplot, p1, width = 8, height = 15)
