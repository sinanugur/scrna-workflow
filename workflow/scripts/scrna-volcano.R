#!/usr/bin/env Rscript


option_list <- list(
    optparse::make_option(c("--rds"),
        type = "character", default = NULL,
        help = "RDS file of marker data frame", metavar = "character"
    ),
    optparse::make_option(c("--vplot"),
        type = "character", default = "vplot.pdf",
        help = "Output volcano plot name", metavar = "character"
    )
)



opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

if (is.null(opt$rds)) {
    optparse::print_help(opt_parser)
    stop("At least one argument must be supplied (rds file)", call. = FALSE)
}

require(EnhancedVolcano)
require(tidyverse)
# try({source("workflow/scripts/scrna-functions.R")})
# try({source(paste0(system("python -c 'import os; import cellsnake; print(os.path.dirname(cellsnake.__file__))'", intern = TRUE),"/scrna/workflow/scripts/scrna-functions.R"))})




all_markers <- readRDS(file = opt$rds)

pdf(opt$vplot, width = 7, height = 7)
for (i in unique(all_markers$condition)) {
    markers <- all_markers %>% filter(condition == i)

    print(EnhancedVolcano(markers, x = "avg_log2FC", y = "p_val_adj", lab = markers %>% pull(gene), FCcutoff = 0.75, subtitle = i))
}
dev.off()