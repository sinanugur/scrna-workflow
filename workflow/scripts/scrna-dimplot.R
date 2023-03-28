#!/usr/bin/env Rscript


option_list <- list(
      optparse::make_option(c("--rds"),
            type = "character", default = NULL,
            help = "Processed rds file of a Seurat object", metavar = "character"
      ),
      optparse::make_option(c("--reduction.type"),
            type = "character", default = "umap",
            help = "Reduction type, umap or tsne", metavar = "character"
      ),
      optparse::make_option(c("--pdfplot"),
            type = "character", default = "reduction.pdf",
            help = "Plot file name", metavar = "character"
      ),
      optparse::make_option(c("--htmlplot"),
            type = "character", default = "reduction.html",
            help = "Plot file name", metavar = "character"
      ),
      optparse::make_option(c("--csv"),
            type = "character", default = NULL,
            help = "CSV meta file", metavar = "character"
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

require(plotly)
require(Seurat)
require(tidyverse)

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


if (!is.null(opt$csv)) {
      metadata <- read.csv(
            opt$csv,
            row.names = 1
      )

      scrna@meta.data <- scrna@meta.data %>%
            tibble::rownames_to_column("barcodes") %>%
            dplyr::left_join(metadata %>% as.data.frame() %>% rownames_to_column("barcodes"), by = "barcodes") %>%
            tibble::column_to_rownames("barcodes")
}

Idents(object = scrna) <- scrna@meta.data[[opt$idents]]

n <- length(Idents(scrna) %>% unique())

palette <- function_color_palette(n)

p1 <- DimPlot(scrna, reduction = opt$reduction.type) & scale_color_manual(values = palette)

ggplotly(p1) -> p1_plotly

p1_plotly %>% htmlwidgets::saveWidget(file = opt$htmlplot, selfcontained = T)



p1 <- DimPlot(scrna, reduction = opt$reduction.type) & scale_color_manual(values = palette)
ggsave(plot = p1, filename = opt$pdfplot, width = 9 + floor(n / 11), height = 7)