#!/usr/bin/env Rscript


option_list = list(
  optparse::make_option(c("--resolution"), type="double", default=0.8, 
              help="Resolution [default= %default]", metavar="character"),
    optparse::make_option(c("--rds"), type="character", default=NULL, 
              help="Processed rds file of a Seurat object", metavar="character"),
    optparse::make_option(c("--logfc.threshold"), type="double", default=0.25, 
              help="LogFC [default= %default]", metavar="character"),
    optparse::make_option(c("--test.use"), type="character", default="wilcox", 
              help="Test use [default= %default]", metavar="character"),
        optparse::make_option(c("--output.xlsx.positive"), type="character", default=NULL, 
              help="Excel table of positive markers", metavar="character"),
    optparse::make_option(c("--output.xlsx.all"), type="character", default=NULL, 
              help="Excel table of all markers", metavar="character")


)
 


opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$rds)){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file)", call.=FALSE)
}

require(Seurat)
require(tidyverse)



source("workflow/scripts/scrna-functions.R")

scrna=readRDS(file = opt$rds)



RNA_=paste0("RNA_snn_res.",opt$resolution)


Idents(object = scrna) <- scrna@meta.data[[RNA_]]


all_markers=FindAllMarkers(scrna, logfc.threshold = opt$logfc.threshold,test.use = opt$test.use )



openxlsx::write.xlsx(all_markers,file=opt$output.xlsx.all)

openxlsx::write.xlsx(all_markers %>% filter(avg_log2FC > 0),file=opt$output.xlsx.positive)
