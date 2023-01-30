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
  optparse::make_option(c("--h5seurat"),
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

if (is.null(opt$h5seurat)) {
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (data.dir and sampleid)", call. = FALSE)
}

require(optparse)
require(tidyverse)
require(Seurat)



# nFeature_RNA is the number of genes detected in each cell. nCount_RNA is the total number of molecules detected within a cell.


#CreateSeuratObject(scrna.data,min.cells = 1,min.features = 5) -> scrna

CreateSeuratObject(LoadH5Seurat(opt$h5seurat)[["RNA"]]@counts,min.cells = opt$min.cells,min.features = opt$min.features)[["RNA"]]@counts %>% as.matrix() %>% 
t() %>% as.data.frame() %>% select(-starts_with("Homo")) %>% saveRDS(.x, file = opt$output.rds)
