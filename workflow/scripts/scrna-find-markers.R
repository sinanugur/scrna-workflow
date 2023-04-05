#!/usr/bin/env Rscript


option_list <- list(
  optparse::make_option(c("--rds"),
    type = "character", default = NULL,
    help = "Processed rds file of a Seurat object", metavar = "character"
  ),
  optparse::make_option(c("--logfc.threshold"),
    type = "double", default = 0.25,
    help = "LogFC [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--test.use"),
    type = "character", default = "wilcox",
    help = "Test use [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--output.rds"),
    type = "character", default = NULL,
    help = "RDS table of all markers", metavar = "character"
  ),
  optparse::make_option(c("--idents"),
    type = "character", default = "seurat_clusters",
    help = "Meta data column name for marker analysis", metavar = "character"
  )
)



opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

if (is.null(opt$rds)) {
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file)", call. = FALSE)
}

require(Seurat)
require(tidyverse)
# try({source("workflow/scripts/scrna-functions.R")})
# try({source(paste0(system("python -c 'import os; import cellsnake; print(os.path.dirname(cellsnake.__file__))'", intern = TRUE),"/scrna/workflow/scripts/scrna-functions.R"))})


scrna <- readRDS(file = opt$rds)
DefaultAssay(scrna) <- "RNA"



# RNA_=paste0("RNA_snn_res.",opt$resolution)
# Idents(object = scrna) <- scrna@meta.data[[RNA_]]
Idents(object = scrna) <- scrna@meta.data[[opt$idents]]


all_markers <- FindAllMarkers(scrna, logfc.threshold = opt$logfc.threshold, test.use = opt$test.use)



saveRDS(all_markers, file = opt$output.rds)