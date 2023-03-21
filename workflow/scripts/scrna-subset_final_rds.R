#!/usr/bin/env Rscript


option_list <- list(
    optparse::make_option(c("--rds"),
        type = "character", default = NULL,
        help = "Processed rds file of a Seurat object", metavar = "character"
    ),
    optparse::make_option(c("--keywords"),
        type = "character", default = NULL,
        help = "Keywords to slice object", metavar = "character"
    ),
    optparse::make_option(c("--idents"),
        type = "character", default = NULL,
        help = "Meta data column name for selecting the slice", metavar = "character"
    ),
    optparse::make_option(c("--output.rds"),
        type = "character", default = "output.rds",
        help = "Output RDS file name [default= %default]", metavar = "character"
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



scrna <- readRDS(file = opt$rds)


if (!is.null(opt$keywords)) {
    keywords <- strsplit(opt$keywords, ",")[[1]]
    patterns <- paste0("(?i)", "(", paste(keywords, collapse = "|"), ")")
    scrna@meta.data %>% dplyr::filter(if_any(everything(), str_detect, pattern = patterns))
}