#!/usr/bin/env Rscript


option_list <- list(
    optparse::make_option(c("--rds"),
        type = "character", default = NULL,
        help = "RDS file of marker data frame", metavar = "character"
    ),
    optparse::make_option(c("--output.xlsx.positive"),
        type = "character", default = NULL,
        help = "Excel table of positive markers", metavar = "character"
    ),
    optparse::make_option(c("--output.xlsx.all"),
        type = "character", default = NULL,
        help = "Excel table of all markers", metavar = "character"
    )
)



opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

if (is.null(opt$rds)) {
    optparse::print_help(opt_parser)
    stop("At least one argument must be supplied (rds file)", call. = FALSE)
}


require(tidyverse)
# try({source("workflow/scripts/scrna-functions.R")})
# try({source(paste0(system("python -c 'import os; import cellsnake; print(os.path.dirname(cellsnake.__file__))'", intern = TRUE),"/scrna/workflow/scripts/scrna-functions.R"))})


all_markers <- readRDS(file = opt$rds)



openxlsx::write.xlsx(all_markers, file = opt$output.xlsx.all)

openxlsx::write.xlsx(all_markers %>% filter(avg_log2FC > 0), file = opt$output.xlsx.positive)