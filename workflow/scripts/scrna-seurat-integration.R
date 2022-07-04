#!/usr/bin/env Rscript

option_list = list(
    optparse::make_option(c("--sampleid"), type="character", default=NULL, 
              help="Sample ID", metavar="character"),
    optparse::make_option(c("-r","--rds"), type="character", default=NULL, 
              help="A list of RDS files of Seurat objects", metavar="character")


)
 
opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$rds) || is.null(opt$sampleid) ){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds and sampleid)", call.=FALSE)
}

require(tidyverse)
require(Seurat)
require(patchwork)
source("workflow/scripts/scrna-functions.R")

files= unlist(strsplit(opt$rds, " "))
print(files)
for(i in files) {

    if(!exists("scrna_list")) {

        scrna_list=list(readRDS(file = i))

    } else {

        scrna_list=append(scrna_list,readRDS(file = i))
    }
    


}



scrna_anchors    <- FindIntegrationAnchors(object.list = scrna_list, dims = 1:30)


scrna     <- IntegrateData(anchorset = scrna_anchors, dims = 1:30)


scrna <- ScaleData(scrna)
scrna <- RunPCA(scrna)
dimensionReduction=function_pca_dimensions(scrna)

scrna <- RunUMAP(scrna, dims = 1:dimensionReduction)

scrna <- FindNeighbors(scrna, reduction = "pca", dims = 1:dimensionReduction)
scrna <- FindClusters(scrna, resolution = c(0.8,2.5))


output.dir=paste0("results/integration/seurat/technicals/")
dir.create(output.dir,recursive = T)

DimPlot(scrna,reduction = "umap",group.by="orig.ident")
ggsave(file=paste0(output.dir,opt$sampleid,"-after-integration-umap.pdf"))




output.dir=paste0("analyses/seuratintegration/")
dir.create(output.dir,recursive = T)

saveRDS(scrna,file = paste0(output.dir,opt$sampleid,"_seurat",".rds"))