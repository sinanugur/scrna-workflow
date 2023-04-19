#!/usr/bin/env Rscript

option_list <- list(
    optparse::make_option(c("--sampleid"),
        type = "character", default = NULL,
        help = "Sample ID", metavar = "character"
    ),
    optparse::make_option(c("-r", "--rds"),
        type = "character", default = NULL,
        help = "A list of RDS files of Seurat objects", metavar = "character"
    ),
    optparse::make_option(c("--output.rds"),
        type = "character", default = "output.rds",
        help = "Output RDS file name [default= %default]", metavar = "character"
    ),
    optparse::make_option(c("--cca.dims"),
        type = "integer", default = 30,
        help = "Which dimensions to use from the CCA to specify the neighbor search space 1 to [default= %default]", metavar = "character"
    ),
    optparse::make_option(c("--reduction"),
        type = "character", default = "cca",
        help = "Integration reduction type [default= %default]", metavar = "character"
    )
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

if (is.null(opt$rds) || is.null(opt$sampleid)) {
    optparse::print_help(opt_parser)
    stop("At least one argument must be supplied (rds and sampleid)", call. = FALSE)
}

require(tidyverse)
require(Seurat)
require(patchwork)
try(
    {
        source("workflow/scripts/scrna-functions.R")
    },
    silent = TRUE
)
try(
    {
        source(paste0(system("python -c 'import os; import cellsnake; print(os.path.dirname(cellsnake.__file__))'", intern = TRUE), "/scrna/workflow/scripts/scrna-functions.R"))
    },
    silent = TRUE
)


files <- unlist(strsplit(opt$rds, " "))
print(files)
for (i in files) {
    if (!exists("scrna_list")) {
        scrna_list <- list(readRDS(file = i))
    } else {
        scrna_list <- append(scrna_list, readRDS(file = i))
    }
}



scrna_anchors <- FindIntegrationAnchors(object.list = scrna_list, dims = 1:opt$cca.dims, reduction = opt$reduction)


scrna <- IntegrateData(anchorset = scrna_anchors, dims = 1:opt$cca.dims)


saveRDS(scrna, file = opt$output.rds)