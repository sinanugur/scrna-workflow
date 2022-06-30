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
require(harmony)
source("workflow/scripts/scrna-functions.R")

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

scrna <- FindVariableFeatures(scrna, selection.method = "vst", nfeatures = 2000)

scrna <- ScaleData(scrna)

scrna <- RunPCA(scrna)

dimensionReduction=function_pca_dimensions(scrna)

scrna <- RunUMAP(scrna, dims = 1:dimensionReduction)



scrna <- scrna %>% RunHarmony("orig.ident",plot_convergence=T)


scrna <- scrna %>% 
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = F) %>% 
  FindNeighbors(reduction = "harmony", k.param = 10, dims = 1:30) %>% 
  FindClusters() %>% 
  identity()

output.dir=paste0("analyses/harmony/")
dir.create(output.dir,recursive = T)

saveRDS(scrna,file = paste0(output.dir,opt$sampleid,"_harmony",".rds"))






