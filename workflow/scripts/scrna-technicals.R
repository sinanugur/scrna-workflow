#!/usr/bin/env Rscript
option_list <- list(
  optparse::make_option(c("--rds"),
    type = "character", default = NULL,
    help = "Processed rds file of a Seurat object", metavar = "character"
  ),
  optparse::make_option(c("--sampleid"),
    type = "character", default = NULL,
    help = "Sample ID", metavar = "character"
  ),
  optparse::make_option(c("--fplot"),
    type = "character", default = NULL,
    help = "nFeature plot", metavar = "character"
  ),
  optparse::make_option(c("--cplot"),
    type = "character", default = NULL,
    help = "nCount plot", metavar = "character"
  ),
  optparse::make_option(c("--mtplot"),
    type = "character", default = NULL,
    help = "Percent MT plot", metavar = "character"
  ),
  optparse::make_option(c("--rpplot"),
    type = "character", default = NULL,
    help = "Ribo plot", metavar = "character"
  )
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

if (is.null(opt$rds)) {
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file and sampleid)", call. = FALSE)
}

require(Seurat)
require(tidyverse)
# try({source("workflow/scripts/scrna-functions.R")})
# try({source(paste0(system("python -c 'import os; import cellsnake; print(os.path.dirname(cellsnake.__file__))'", intern = TRUE),"/scrna/workflow/scripts/scrna-functions.R"))})



scrna <- readRDS(file = opt$rds)


FeaturePlot(scrna, features = "nFeature_RNA", pt.size = 0.1)

ggsave(opt$fplot, width = 7, height = 5)



FeaturePlot(scrna, features = "nCount_RNA", pt.size = 0.1)

ggsave(opt$cplot, width = 7, height = 5)


FeaturePlot(scrna, features = "percent.mt", pt.size = 0.1)

ggsave(opt$mtplot, width = 7, height = 5)

FeaturePlot(scrna, features = "percent.rp", pt.size = 0.1)

ggsave(opt$rpplot, width = 7, height = 5)