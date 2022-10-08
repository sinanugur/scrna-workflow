#!/usr/bin/env Rscript


option_list = list(
  optparse::make_option(c("--resolution"), type="double", default=0.8, 
              help="Resolution [default= %default]", metavar="character"),

    optparse::make_option(c("--rds"), type="character", default=NULL, 
              help="Processed rds file of a Seurat object", metavar="character"),
      optparse::make_option(c("--output.umap.plot"), type="character", default="umap.pdf", 
              help="UMAP plot file name", metavar="character")


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

p1 <- DimPlot(scrna, reduction = "umap", label = TRUE,label.size = 10) 





ggsave(plot =p1,filename=opt$output.umap.plot,width=13,height=7)

