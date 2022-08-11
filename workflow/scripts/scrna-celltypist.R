#!/usr/bin/env Rscript

option_list = list(
    optparse::make_option(c("--sampleid"), type="character", default=NULL, 
              help="Sample ID", metavar="character"),
    optparse::make_option(c("-r","--rds"), type="character", default=NULL, 
              help="A list of RDS files of Seurat objects", metavar="character"),
        optparse::make_option(c("--csv"), type="character", default=NULL, 
              help="Celltypist prediction file", metavar="character"),
        optparse::make_option(c("-o","--output"), type="character", default=NULL, 
              help="Output RDS file", metavar="character")

)
 
opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$rds) || is.null(opt$sampleid) || is.null(opt$csv) || is.null(opt$output) ){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds and sampleid)", call.=FALSE)
}

require(tidyverse)
require(Seurat)
source("workflow/scripts/scrna-functions.R")


scrna=readRDS(file = opt$rds)

celltyp=read.csv(
    opt$csv,row.names = 1)

scrna@meta.data <- scrna@meta.data %>% tibble::rownames_to_column("barcodes") %>% dplyr::left_join(celltyp %>% as.data.frame() %>% rownames_to_column("barcodes"), by = "barcodes") %>% tibble::column_to_rownames("barcodes")


saveRDS(scrna,file = opt$output)