#!/usr/bin/env Rscript

option_list <- list(
      optparse::make_option(c("--rds"),
            type = "character", default = NULL,
            help = "A list of RDS files of Seurat objects", metavar = "character"
      ),
      optparse::make_option(c("--sheplot"),
            type = "character", default = "sheplot.pdf",
            help = "Output score heatmap plot file name", metavar = "character"
      ),
      optparse::make_option(c("--sheplottop"),
            type = "character", default = "sheplot.pdf",
            help = "Output score heatmap plot file name, top 20", metavar = "character"
      ),
      optparse::make_option(c("--pheplot"),
            type = "character", default = "pheplot.pdf",
            help = "Output heatmap plot file name", metavar = "character"
      ),
      optparse::make_option(c("--idents"),
            type = "character", default = "seurat_clusters",
            help = "Meta data column name", metavar = "character"
      ),
      optparse::make_option(c("--csv"),
            type = "character", default = NULL,
            help = "A meta data table", metavar = "character"
      ),
      optparse::make_option(c("--prediction"),
            type = "character", default = "pred.rds",
            help = "Input prediction file", metavar = "character"
      ),
)


opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

if (is.null(opt$rds)) {
      optparse::print_help(opt_parser)
      stop("At least one argument must be supplied (rds and sampleid)", call. = FALSE)
}


require(optparse)
require(SingleR)
# require(celldex)
require(pheatmap)
require(Seurat)

scrna <- readRDS(file = opt$rds)
DefaultAssay(scrna) <- "RNA"



# celltype annotation with SingleR
pred <- readRDS(opt$prediction)


plotScoreHeatmap(pred) -> p1
ggsave(plot = p1, filename = opt$sheplot, width = 15, height = 8)


plotScoreHeatmap(pred, show.labels = F, max.labels = 20) -> p1
ggsave(plot = p1, filename = opt$sheplottop, width = 7, height = 4)


tab <- table(Assigned = pred$pruned.labels, Cluster = scrna[[opt$idents]])


pheatmap(log2(tab + 10), color = colorRampPalette(c("white", "blue"))(101)) -> p1



n1 <- length(unique(scrna$seurat_clusters))
n2 <- length(unique(pred$pruned.labels))

ggsave(plot = p1, filename = opt$pheplot, width = 6 + (n1 * 0.10), height = 4 + (n2 * 0.10))