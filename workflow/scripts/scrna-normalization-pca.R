#!/usr/bin/env Rscript


option_list = list(
  optparse::make_option(c("--scale.factor"), type="integer", default=10000, 
              help="Scale factor [default= %default]", metavar="integer"),
  optparse::make_option(c("--nfeatures"), type="integer", default=2000, 
              help="Highly variable features [default= %default]", metavar="integer"),
    optparse::make_option(c("--rds"), type="character", default=NULL, 
              help="RAW rds file of a Seurat object", metavar="character"),
    optparse::make_option(c("--sampleid"), type="character", default=NULL, 
              help="Sample ID", metavar="character"),
    optparse::make_option(c("--normalization.method"), type="character", default="LogNormalize", 
              help="Normalization method[default= %default]", metavar="character"),
    optparse::make_option(c("--minCov"), type="integer", default=1000, 
              help="if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount [default= %default]", metavar="character")


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


scrna=readRDS(file = opt$rds)


scrna <- NormalizeData(scrna, normalization.method = opt$normalization.method, scale.factor = opt$scale.factor)

scrna <- FindVariableFeatures(scrna, selection.method = "vst", nfeatures = 2000)

output.dir=paste0("results/",opt$sampleid,"/technicals/")
dir.create(output.dir,recursive = T)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(scrna), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(scrna)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)



ggsave(paste0(output.dir,"highly-variable-features.pdf"), plot2 ,width = 8,height = 9)