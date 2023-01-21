#!/usr/bin/env Rscript

option_list <- list(
  optparse::make_option(c("--min.cells"),
    type = "integer", default = 1,
    help = "Min cells [default= %default]", metavar = "integer"
  ),
  optparse::make_option(c("--min.features"),
    type = "integer", default = 5,
    help = "Min features [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--data.dir"),
    type = "character", default = NULL,
    help = "Data directory", metavar = "character"
  ),
  optparse::make_option(c("--sampleid"),
    type = "character", default = NULL,
    help = "Sample ID", metavar = "character"
  ),
  optparse::make_option(c("--output.rds"),
    type = "character", default = "output.rds",
    help = "Output RDS file name [default= %default]", metavar = "character"
  )

)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

if (is.null(opt$data.dir) || is.null(opt$sampleid)) {
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (data.dir and sampleid)", call. = FALSE)
}

require(optparse)
require(tidyverse)
require(Seurat)
require(patchwork)
require(tools)
require(data.table)


# nFeature_RNA is the number of genes detected in each cell. nCount_RNA is the total number of molecules detected within a cell.




try({
  scrna.data <- Read10X(data.dir = opt$data.dir)
},silent = TRUE)
try({
  scrna.data <- Read10X_h5(filename = paste0(opt$data.dir, "/filtered_feature_bc_matrix.h5"))
},silent = TRUE)
try({
  scrna.data <- Read10X_h5(filename = opt$data.dir)
},silent = TRUE)




CreateSeuratObject(scrna.data,min.cells = 1,min.features = 5) -> scrna


saveRDS(scrna, file = opt$output.rds)
