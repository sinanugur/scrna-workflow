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
  optparse::make_option(c("--idents"),
    type = "character", default = "seurat_clusters",
    help = "Meta data column name for marker analysis", metavar = "character"
  ),
  optparse::make_option(c("--ccplot"),
    type = "character", default = "ccplot.pdf",
    help = "Cell cluster count plot", metavar = "character"
  ),
  optparse::make_option(c("--ccbarplot"),
    type = "character", default = "ccbarplot.pdf",
    help = "Cell cluster count plot", metavar = "character"
  ),
  optparse::make_option(c("--htmlplot"),
    type = "character", default = "htmlplot.pdf",
    help = "Cell cluster html plot", metavar = "character"
  )
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

if (is.null(opt$rds)) {
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file and sampleid)", call. = FALSE)
}

require(plotly)
require(ggpubr)
require(Seurat)
require(tidyverse)
# try({source("workflow/scripts/scrna-functions.R")})
# try({source(paste0(system("python -c 'import os; import cellsnake; print(os.path.dirname(cellsnake.__file__))'", intern = TRUE),"/scrna/workflow/scripts/scrna-functions.R"))})



scrna <- readRDS(file = opt$rds)



scrna@meta.data %>%
  dplyr::add_count(orig.ident) %>%
  dplyr::mutate(total_clusters = length(unique(get(opt$idents)))) %>%
  distinct(orig.ident, n, total_clusters) %>%
  dplyr::select("Sample Name" = orig.ident, "Total Cells" = n, "Total Clusters" = total_clusters) %>%
  ggtexttable(rows = NULL, theme = ttheme("light")) -> p1

ggsave(opt$ccplot, p1)

scrna@meta.data %>%
  dplyr::group_by(orig.ident, get(opt$idents)) %>%
  dplyr::count() %>%
  dplyr::select("Sample Name" = 1, "Cluster" = 2, "Total Cells" = 3) %>%
  ggplot(aes(x = Cluster, y = `Total Cells`, fill = `Sample Name`)) +
  geom_col() +
  ggthemes::theme_hc() +
  theme(legend.title = element_blank()) -> p2
n <- length(unique(scrna@meta.data[opt$idents]))

ggsave(opt$ccbarplot, p2, height = 5, width = 5 + (n * 0.12))


ggplotly(p2) -> p1_plotly

p1_plotly %>% htmlwidgets::saveWidget(file = opt$htmlplot, selfcontained = T)