#!/usr/bin/env Rscript

option_list <- list(
  optparse::make_option(c("--rds"),
    type = "character", default = NULL,
    help = "A list of RDS files of Seurat objects", metavar = "character"
  ),
  optparse::make_option(c("--output"),
    type = "character", default = "output.h5ad",
    help = "Output h5ad file name", metavar = "character"
  )
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

if (is.null(opt$rds)) {
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file and sampleid)", call. = FALSE)
}

if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
  remotes::install_github("mojaveazure/seurat-disk")
}
require(Seurat)
require(SeuratDisk)
require(tidyverse)


scrna <- readRDS(file = opt$rds)

UpdateSeuratObject(scrna) -> scrna
DefaultAssay(scrna) <- "RNA"

DietSeurat(scrna) -> scrna

scrna@meta.data %>% dplyr::mutate(dplyr::across(where(is.factor), as.character)) -> scrna@meta.data


output_file_name <- str_remove_all(opt$output, ".h5ad$")

SaveH5Seurat(scrna, filename = paste0(output_file_name, ".h5Seurat"), overwrite = TRUE)
SeuratDisk::Convert(paste0(output_file_name, ".h5Seurat"), dest = "h5ad", overwrite = TRUE)