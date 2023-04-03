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

if (!requireNamespace("SeuratWrappers", quietly = TRUE)) {
    remotes::install_github("satijalab/seurat-wrappers")
}

if (!requireNamespace("monocle3", quietly = TRUE)) {
    remotes::install_github("cole-trapnell-lab/monocle3")
}

require(monocle3)
require(Seurat)
require(SeuratWrappers)
require(tidyverse)
require(patchwork)
require(magrittr)


if (is.null(opt$rds)) {
    optparse::print_help(opt_parser)
    stop("At least one argument must be supplied (rds file)", call. = FALSE)
}
scrna <- readRDS(file = opt$rds)


cds <- as.cell_data_set(scrna)
cds <- cluster_cells(cds)

p1 <- plot_cells(cds, color_cells_by = "singler", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)


wrap_plots(p1, p2)

ggsave(opt$pplot, width = 11, height = 5.5)


as.Seurat(cds, assay = NULL) -> scrna
scrna@meta.data %>%
    dplyr::count(monocle3_partitions) %>%
    dplyr::mutate(perc = (n * 100) / sum(n)) %>%
    dplyr::filter(perc >= 5, n > 200) %>%
    pull(monocle3_partitions) %>%
    as.vector() -> partitions

print(partitions)
for (i in partitions) {
    integrated.sub <- subset(scrna, monocle3_partitions == i)


    cds2 <- as.cell_data_set(integrated.sub)
    if (i != "1") {
        cds2 <- cluster_cells(cds2)
    }
    cds2 <- learn_graph(cds2)
    p1 <- plot_cells(cds2, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

    # max.avp <- which.max(unlist(FetchData(integrated.sub, "AVP")))
    # max.avp <- colnames(integrated.sub)[max.avp]
    # cds2 <- order_cells(cds2, root_cells = max.avp)
    # p2 <- plot_cells(cds2, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,label_branch_points = FALSE)

    # wrap_plots(p1, p2)
    ggsave(paste0(opt$output.dir, "/plot_monocle-partition-", i, ".pdf"), plot = p1, width = 6, height = 5.5)
}