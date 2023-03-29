#!/usr/bin/env Rscript


# !/usr/bin/env Rscript
option_list <- list(
    optparse::make_option(c("--rds"),
        type = "character", default = NULL,
        help = "Processed rds file of a Seurat object", metavar = "character"
    ),
    optparse::make_option(c("--xlsx"),
        type = "character", default = NULL,
        help = "Excel table of markers for input", metavar = "character"
    ),
    optparse::make_option(c("--output.plot"),
        type = "character", default = "output.pdf",
        help = "Output plot directory", metavar = "character"
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

Idents(object = scrna) <- scrna@meta.data[[opt$idents]]

n <- length(Idents(scrna) %>% unique())

Positive_Features <- openxlsx::read.xlsx(opt$xlsx)


Positive_Features %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC) -> df


df %>%
    distinct(gene) %>%
    pull() -> not.all.genes

scrna <- ScaleData(scrna, features = not.all.genes)

DoHeatmap(object = scrna, features = not.all.genes, label = F) & theme(axis.text.y = element_text(size = 5)) & scale_fill_continuous(type = "viridis") -> p1

ggsave(opt$output.plot, p1, height = 4 + (n * 0.2), width = 5 + (n * 0.05), useDingbats = TRUE)