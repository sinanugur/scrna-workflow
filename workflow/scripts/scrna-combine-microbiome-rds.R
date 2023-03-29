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


files <- unlist(strsplit(opt$rds, " "))
print(files)
for (i in files) {
    if (!exists("scrna")) {
        scrna <- readRDS(file = i)
    } else {
        scrna <- bind_rows(scrna, readRDS(file = i))
    }
}


saveRDS(scrna, file = opt$output.rds)