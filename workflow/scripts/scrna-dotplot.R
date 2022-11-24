#!/usr/bin/env Rscript


option_list = list(
  optparse::make_option(c("--resolution"), type="double", default=0.8, 
              help="Resolution [default= %default]", metavar="character"),
    optparse::make_option(c("--rds"), type="character", default=NULL, 
              help="Processed rds file of a Seurat object", metavar="character"),
        optparse::make_option(c("--reduction.type"), type="character", default="umap", 
              help="Reduction type, umap or tsne", metavar="character"),
        optparse::make_option(c("--tsv"), type="character", default=NULL, 
              help="A text file contains the gene list", metavar="character"),
  optparse::make_option(c("--output.dotplot"), type="character", default=NULL, 
              help="Output dotplot name", metavar="character")



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

markers=read_tsv(opt$tsv,col_names=FALSE) %>% pull()

markers=sort(intersect(markers,rownames(scrna)))

DotPlot(scrna, features = markers, dot.scale = 8) + RotatedAxis()
ggsave(opt$output.dotplot,width = 13,height = 8)

