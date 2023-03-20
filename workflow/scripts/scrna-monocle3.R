#!/usr/bin/env Rscript
option_list <- list(
    optparse::make_option(c("--rds"),
        type = "character", default = NULL,
        help = "Processed rds file of a Seurat object", metavar = "character"
    ),
    optparse::make_option(c("--pplot"),
        type = "character", default = NULL,
        help = "Partition plot", metavar = "character"
    ),
    optparse::make_option(c("--output.dir"),
        type = "character", default = NULL,
        help = "Output plot directory", metavar = "character"
    )
)



opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)


library(monocle3)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(tidyverse)
library(patchwork)
library(magrittr)


if (is.null(opt$rds)) {
    optparse::print_help(opt_parser)
    stop("At least one argument must be supplied (rds file)", call. = FALSE)
}
scrna <- readRDS(file = opt$rds)


cds <- as.cell_data_set(scrna)
cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

ggsave(opt$pplot, width = 10, height = 5)


as.Seurat(cds, assay = NULL) -> scrna
scrna@meta.data %>%
    dplyr::count(monocle3_partitions) %>%
    dplyr::mutate(perc = n * 100 / sum(n)) %>%
    dplyr::filter(perc >= 5) %>%
    pull(monocle3_partitions) %>%
    as.vector() -> partitions


for (i in partitions) {
    integrated.sub <- subset(as.Seurat(cds, assay = NULL), monocle3_partitions == i)

    cds <- as.cell_data_set(integrated.sub)
    cds <- learn_graph(cds)
    plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

    max.avp <- which.max(unlist(FetchData(integrated.sub, "AVP")))
    max.avp <- colnames(integrated.sub)[max.avp]
    cds <- order_cells(cds, root_cells = max.avp)
    plot_cells(cds,
        color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,
        label_branch_points = FALSE
    )
    ggsave(paste0(opt$output.plot.dir, "/plot_monocle3-partition-", i, ".pdf"), width = 6, height = 6)
}