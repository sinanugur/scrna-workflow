#!/usr/bin/env Rscript

option_list = list(
    optparse::make_option(c("--sampleid"), type="character", default=NULL, 
              help="Sample ID", metavar="character"),
    optparse::make_option(c("-r","--rds"), type="character", default=NULL, 
              help="A list of RDS files of Seurat objects", metavar="character"),
    optparse::make_option(c("--resolution"), type="double", default=0.8, 
              help="Resolution [default= %default]", metavar="character")


)

opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$rds) || is.null(opt$sampleid) ){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file and sampleid)", call.=FALSE)
}

require(Seurat)
require(SeuratDisk)
require(tidyverse)


scrna=readRDS(file = opt$rds)

UpdateSeuratObject(scrna) -> scrna
DefaultAssay(scrna) <- "RNA"

DietSeurat(scrna) -> scrna

scrna@meta.data %>% dplyr::mutate(dplyr::across(where(is.factor), as.character)) -> scrna@meta.data



output.dir=paste0("analyses/h5ad/",opt$resolution,"/")
dir.create(output.dir,recursive = T)

SaveH5Seurat(scrna,paste0(output.dir,opt$sampleid,".h5Seurat"))
SeuratDisk::Convert(paste0(output.dir,opt$sampleid,".h5Seurat"), dest = "h5ad")