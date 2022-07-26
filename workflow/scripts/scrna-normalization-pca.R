#!/usr/bin/env Rscript


option_list = list(
  optparse::make_option(c("--scale.factor"), type="integer", default=10000, 
              help="Scale factor [default= %default]", metavar="character"),
  optparse::make_option(c("--nfeatures"), type="integer", default=2000, 
              help="Highly variable features [default= %default]", metavar="integer"),
    optparse::make_option(c("--rds"), type="character", default=NULL, 
              help="RAW rds file of a Seurat object", metavar="character"),
    optparse::make_option(c("--sampleid"), type="character", default=NULL, 
              help="Sample ID", metavar="character"),
    optparse::make_option(c("--normalization.method"), type="character", default="LogNormalize", 
              help="Normalization method[default= %default]", metavar="character")


)
 
opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$rds) || is.null(opt$sampleid) ){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file and sampleid)", call.=FALSE)
}

require(optparse)
require(tidyverse)
require(Seurat)
require(patchwork)
require(DoubletFinder)

source("workflow/scripts/scrna-functions.R")

scrna=readRDS(file = opt$rds)


scrna <- NormalizeData(scrna, normalization.method = opt$normalization.method, scale.factor = opt$scale.factor)

scrna <- FindVariableFeatures(scrna, selection.method = "vst", nfeatures = opt$nfeatures)

output.dir=paste0("results/",opt$sampleid,"/technicals/")
dir.create(output.dir,recursive = T)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(scrna), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(scrna)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)



ggsave(paste0(output.dir,"highly-variable-features.pdf"), plot2 ,width = 8,height = 9)



all.genes <- rownames(scrna)
scrna <- ScaleData(scrna, features = all.genes)


scrna <- RunPCA(scrna, features = VariableFeatures(object = scrna))


DimHeatmap(scrna, dims = 1:15, cells = 500, balanced = TRUE,fast = FALSE)

ggsave(paste0(output.dir,"DimHeatMap_plot.pdf") ,width = 8,height = 15)


# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
scrna <- JackStraw(scrna, num.replicate = 100,  dims=50)
scrna <- ScoreJackStraw(scrna, dims = 1:50)


plot1 <- JackStrawPlot(scrna, dims = 1:50) 
plot2 <- ElbowPlot(scrna, ndims=50)

ggsave(paste0(output.dir,"JackandElbow_plot.pdf"), plot1 + plot2,width = 13,height = 5)

dimensionReduction=function_pca_dimensions(scrna)


scrna <- FindNeighbors(scrna, dims = 1:dimensionReduction)
scrna <- FindClusters(scrna, resolution = c(2.5,0.8))

scrna <- RunUMAP(scrna, dims = 1:dimensionReduction)


Idents(object = scrna) <- scrna@meta.data[["RNA_snn_res.0.8"]]

scrna$seurat_clusters <- scrna@meta.data[["RNA_snn_res.0.8"]]



homotypic.prop <- modelHomotypic(scrna$seurat_clusters)         
nExp_poi <- round(0.075*nrow(scrna@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
scrna <- doubletFinder_v3(scrna, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
scrna <- doubletFinder_v3(scrna, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = paste0("pANN_0.25_0.09_",nExp_poi), sct = FALSE)


scrna@meta.data %>% tibble::rownames_to_column("barcodes") %>% select(barcodes,starts_with("DF")) %>% select(barcodes,DoubletFinder=2) -> Doublet_Df

scrna@meta.data <- scrna@meta.data %>% select(!starts_with("DF")) %>% select(!starts_with("pANN")) %>%
  tibble::rownames_to_column("barcodes") %>%
  dplyr::left_join(Doublet_Df, by = "barcodes") %>%
  tibble::column_to_rownames("barcodes")


output.dir=paste0("analyses/processed/")
dir.create(output.dir,recursive = T)

saveRDS(scrna,file = paste0(output.dir,opt$sampleid,".rds"))


