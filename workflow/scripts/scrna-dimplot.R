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
      optparse::make_option(c("--output.reduction.plot"),
            type = "character", default = "reduction.pdf",
            help = "Plot file name", metavar = "character"
      ),
      optparse::make_option(c("--csv"),
            type = "character", default = NULL,
            help = "Celltypist prediction file", metavar = "character"
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
require(randomcoloR)


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
set.seed(149)
palette <- sort(distinctColorPalette(n))


names(palette) <- Idents(scrna) %>%
      unique() %>%
      sort() %>%
      as.character()
print(palette)

p1 <- DimPlot(scrna, reduction = opt$reduction.type, label = TRUE, repel = TRUE) & ggthemes::theme_few() & scale_color_manual(values = palette)




ggsave(plot = p1, filename = opt$output.reduction.plot, width = 9, height = 7)