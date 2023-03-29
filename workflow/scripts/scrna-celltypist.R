#!/usr/bin/env Rscript

option_list <- list(
      optparse::make_option(c("--sampleid"),
            type = "character", default = NULL,
            help = "Sample ID", metavar = "character"
      ),
      optparse::make_option(c("--rds"),
            type = "character", default = NULL,
            help = "A list of RDS files of Seurat objects", metavar = "character"
      ),
      optparse::make_option(c("--csv"),
            type = "character", default = NULL,
            help = "Celltypist prediction file", metavar = "character"
      ),
      optparse::make_option(c("--output.tsne.plot"),
            type = "character", default = NULL,
            help = "Output tsne file name", metavar = "character"
      ),
      optparse::make_option(c("--output.umap.plot"),
            type = "character", default = NULL,
            help = "Output umap file name", metavar = "character"
      ),
      optparse::make_option(c("--percentage"),
            type = "double", default = 5,
            help = "Cluster mimnimum percentage to plot", metavar = "double"
      ), optparse::make_option(c("--labels"), action = "store_true", default = FALSE, help = "Print labels on the plot")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

if (is.null(opt$rds) || is.null(opt$csv)) {
      optparse::print_help(opt_parser)
      stop("At least one argument must be supplied (rds and sampleid)", call. = FALSE)
}

require(tidyverse)
require(Seurat)


tryCatch(
      {
            source("workflow/scripts/scrna-functions.R")
      },
      error = function(cond) {
            source(paste0(system("python -c 'import os; import cellsnake; print(os.path.dirname(cellsnake.__file__))'", intern = TRUE), "/scrna/workflow/scripts/scrna-functions.R"))
      }
)


scrna <- readRDS(file = opt$rds)
DefaultAssay(scrna) <- "RNA"

celltypist <- read.csv(
      opt$csv,
      row.names = 1
)

scrna@meta.data <- scrna@meta.data %>%
      tibble::rownames_to_column("barcodes") %>%
      dplyr::left_join(celltypist %>% as.data.frame() %>% rownames_to_column("barcodes"), by = "barcodes") %>%
      tibble::column_to_rownames("barcodes")


n <- length(scrna@meta.data %>% pull(majority_voting) %>% unique())

palette <- function_color_palette(n)
palette <- function_color_palette(n)
palette <- setNames(palette, scrna@meta.data %>% pull(majority_voting) %>% unique())

breaks <- scrna@meta.data %>%
      dplyr::select(majority_voting) %>%
      dplyr::count(majority_voting) %>%
      dplyr::mutate(perc = (n * 100) / sum(n)) %>%
      dplyr::filter(perc >= opt$percentage) %>%
      dplyr::select(-n, -perc) %>%
      distinct() %>%
      pull() %>%
      as.character()

p1 <- DimPlot(scrna, reduction = "tsne", label = opt$labels, group.by = "majority_voting", repel = TRUE) & scale_color_manual(values = palette, breaks = breaks)
p2 <- DimPlot(scrna, reduction = "umap", label = opt$labels, group.by = "majority_voting", repel = TRUE) & scale_color_manual(values = palette, breaks = breaks)


m <- max(str_count(breaks))

w <- c(7.5 + (m * 0.08) * (floor(length(breaks) / 11) + 1))


ggsave(plot = p1, filename = opt$output.tsne.plot, width = w, height = 7)
ggsave(plot = p2, filename = opt$output.umap.plot, width = w, height = 7)



# saveRDS(scrna,file = opt$output)