#!/usr/bin/env Rscript
option_list <- list(
    optparse::make_option(c("--rds"),
        type = "character", default = NULL,
        help = "Processed rds file of a Seurat object", metavar = "character"
    ),
    optparse::make_option(c("--top_n"),
        type = "integer", default = 50,
        help = "How many features to plot per cluster [default= %default]", metavar = "integer"
    ),
    optparse::make_option(c("--xlsx"),
        type = "character", default = NULL,
        help = "Excel table of markers for input", metavar = "character"
    ),
    optparse::make_option(c("--output.plot.dir"),
        type = "character", default = NULL,
        help = "Output plot directory", metavar = "character"
    ),
    optparse::make_option(c("--reduction.type"),
        type = "character", default = "umap",
        help = "Reduction type, umap or tsne", metavar = "character"
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
require(viridis)
require(randomcoloR)




scrna <- readRDS(file = opt$rds)

n<-length(Idents(scrna) %>% unique())
set.seed(149)
palette <- distinctColorPalette(n)


Positive_Features <- openxlsx::read.xlsx(opt$xlsx) %>%
    group_by(cluster) %>%
    slice_min(order_by = p_val_adj, n = opt$top_n) %>%
    select(cluster, gene)


for (d in (Positive_Features %>% distinct(cluster) %>% pull())) {
    dir.create(paste0(opt$output.plot.dir, "/cluster", d, "/"), recursive = TRUE)
}

options(warn = -1)




suppressMessages(for (i in 1:nrow(Positive_Features)) {
    gene <- Positive_Features[i, ]$gene
    cluster <- Positive_Features[i, ]$cluster

    p1 <- FeaturePlot(scrna, reduction = opt$reduction.type, features = gene) & theme_cellsnake_classic() & scale_color_continuous(type = "viridis") & labs(color="Expression") & theme(axis.text = element_text(size=12))
    p2 <- DotPlot(scrna, features = gene) & theme_cellsnake_classic() & scale_color_continuous(type = "viridis") & labs(color="Average Expression",size="Percent Expressed") & ylab("Identity") &  theme(axis.title.x = element_blank(),axis.text = element_text(size=12)) & theme(legend.position = "right")
    p3 <- VlnPlot(scrna, features = gene) & theme_cellsnake_classic() & scale_fill_manual(values = palette) & theme(legend.position = "right",axis.text = element_text(size = 12)) & labs(fill="") & xlab("Identity") & ylab("Expression Level")

    suppressWarnings(((p1 | p2) / p3) -> wp)

    ggsave(paste0(opt$output.plot.dir, "/cluster", cluster, "/", gene, ".pdf"), wp,height=6+(n*0.15),width=7+(n*0.15),useDingbats = TRUE)
})