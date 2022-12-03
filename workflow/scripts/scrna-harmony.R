#!/usr/bin/env Rscript

option_list = list(
    optparse::make_option(c("--sampleid"), type="character", default=NULL, 
              help="Sample ID", metavar="character"),
    optparse::make_option(c("-r","--rds"), type="character", default=NULL, 
              help="A list of RDS files of Seurat objects", metavar="character"),
    optparse::make_option(c("--nfeatures"), type="integer", default=2000, 
              help="Highly variable features [default= %default]", metavar="integer"),
    optparse::make_option(c("--normalization.method"), type="character", default="LogNormalize", 
              help="Normalization method[default= %default]", metavar="character"),
    optparse::make_option(c("--scale.factor"), type="integer", default=10000, 
              help="Scale factor [default= %default]", metavar="character"),
    optparse::make_option(c("--resolution"), type="double", default=0.8, 
              help="Resolution [default= %default]", metavar="character")


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
require(harmony)
try({source("workflow/scripts/scrna-functions.R")})
try({source(paste0(system("python -c 'import os; import cellsnake; print(os.path.dirname(cellsnake.__file__))'", intern = TRUE),"/scrna/workflow/scripts/scrna-functions.R"))})


files= unlist(strsplit(opt$rds, " "))
print(files)
for(i in files) {

    if(!exists("scrna")) {

        scrna=readRDS(file = i)

    } else {

        scrna=merge(scrna,readRDS(file = i))
    }
    


}


scrna <- NormalizeData(scrna, normalization.method = opt$normalization.method, scale.factor = opt$scale.factor)

scrna <- FindVariableFeatures(scrna, selection.method = "vst", nfeatures = opt$nfeatures)

scrna <- ScaleData(scrna)

scrna <- RunPCA(scrna)

dimensionReduction=function_pca_dimensions(scrna)

scrna <- RunUMAP(scrna, dims = 1:dimensionReduction)



output.dir=paste0("results/integration/harmony/technicals/")
dir.create(output.dir,recursive = T)
DimPlot(scrna,reduction = "umap",group.by="orig.ident")
ggsave(file=paste0(output.dir,opt$sampleid,"-before-integration-umap.pdf"))





scrna <- scrna %>% RunHarmony("orig.ident",plot_convergence=T)


scrna <- scrna %>% 
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = F) %>% 
  FindNeighbors(reduction = "harmony", k.param = 10, dims = 1:30) %>% 
  FindClusters(resolution = opt$resolution) %>% 
  identity()


DimPlot(scrna,reduction = "umap",group.by="orig.ident")
ggsave(file=paste0(output.dir,opt$sampleid,"-after-integration-umap.pdf"))



output.dir=paste0("analyses/integration/harmony/")
dir.create(output.dir,recursive = T)

saveRDS(scrna,file = paste0(output.dir,opt$sampleid,"_harmony",".rds"))






