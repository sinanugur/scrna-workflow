#!/usr/bin/env Rscript


option_list = list(
    optparse::make_option(c("--rds"), type="character", default=NULL, 
              help="RAW rds file of a Seurat object", metavar="character"),
    optparse::make_option(c("--sampleid"), type="character", default=NULL, 
              help="Sample ID", metavar="character")


)
 
opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$rds) || is.null(opt$sampleid) ){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file and sampleid)", call.=FALSE)
}

require(optparse)
require(Seurat)
require(clustree)


source("workflow/scripts/scrna-functions.R")

scrna=readRDS(file = opt$rds)


scrna <- FindClusters(scrna, resolution = seq(0.1,2.5,0.1))

clustree(scrna) -> p1

output.dir=paste0("results/",opt$sampleid,"/clusteringTree/")
dir.create(output.dir,recursive = T)

ggsave(paste0(output.dir,"/clusteringTree-",opt$sampleid,".pdf"),p1,width=8,height=15)