#!/usr/bin/env Rscript
option_list = list(
  optparse::make_option(c("--resolution"), type="double", default=0.8, 
              help="Resolution [default= %default]", metavar="character"),

    optparse::make_option(c("--rds"), type="character", default=NULL, 
              help="Processed rds file of a Seurat object", metavar="character"),

    optparse::make_option(c("--sampleid"), type="character", default=NULL, 
              help="Sample ID", metavar="character")


)
 


opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$rds) || is.null(opt$sampleid) ){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file and sampleid)", call.=FALSE)
}

require(Seurat)
require(tidyverse)



source("workflow/scripts/scrna-functions.R")


scrna=readRDS(file = opt$rds)

RNA_=paste0("RNA_snn_res.",opt$resolution)

metrics=table(scrna@meta.data[[RNA_]], scrna@meta.data$orig.ident)

output.dir=paste0("results/",opt$sampleid,"/resolution-",opt$resolution,"/")
dir.create(output.dir,recursive = T)

openxlsx::write.xlsx(metrics %>% as.data.frame() %>% select(Cluster=1,everything()),file=paste0(output.dir,opt$sampleid,".number-of-cells-per-cluster",".xlsx"))