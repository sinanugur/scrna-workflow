#!/usr/bin/env Rscript
option_list <- list(
  optparse::make_option(c("--rds"),
    type = "character", default = NULL,
    help = "Processed rds file of a Seurat object", metavar = "character"
  ),
  optparse::make_option(c("--tsv"),
    type = "character", default = NULL,
    help = "A text file contains the gene list", metavar = "character"
  ),
  optparse::make_option(c("--gene"),
    type = "character", default = NULL,
    help = "A list of genes to create plots", metavar = "character"
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
    help = "Meta data column name", metavar = "character"
  )
)



opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

try(
  {
    source("workflow/scripts/scrna-functions.R")
  },
  silent = TRUE
)
try(
  {
    source(paste0(system("python -c 'import os; import cellsnake; print(os.path.dirname(cellsnake.__file__))'", intern = TRUE), "/scrna/workflow/scripts/scrna-functions.R"))
  },
  silent = TRUE
)

if (is.null(opt$rds)) {
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file)", call. = FALSE)
}

require(Seurat)
require(tidyverse)
require(viridis)

markers <- c()

if (!is.null(opt$tsv)) {
  markers <- read_tsv(opt$tsv, col_names = FALSE) %>% pull()
}



scrna <- readRDS(file = opt$rds)
DefaultAssay(scrna) <- "RNA"

# RNA_=paste0("RNA_snn_res.",opt$resolution)


Idents(object = scrna) <- scrna@meta.data[[opt$idents]]


if (!is.null(opt$gene)) {
  markers <- c(unlist(strsplit(opt$gene, " ")), markers)
}

# identification=opt$idents



# FetchData(scrna,c(identification,"UMAP_1","UMAP_2","TSNE_1","TSNE_2",markers)) %>% rownames_to_column("barcodes") %>% gather(gene,expr,5:last_col()) %>% group_by(gene,across(identification)) %>% dplyr::rename(clusters:=identification) %>% ungroup() -> plot_df


# print(head(plot_df))
dir.create(opt$output.plot.dir, recursive = TRUE)

suppressMessages(for (i in markers) {
  n <- length(Idents(scrna) %>% unique())
  palette <- function_color_palette(n)

  tryCatch(
    {
      p1 <- FeaturePlot(scrna, reduction = opt$reduction.type, features = i) & scale_color_continuous(type = "viridis") & labs(color = "Expression") & theme(axis.text = element_text(size = 10))
      p2 <- DotPlot(scrna, features = i) & scale_color_continuous(type = "viridis") & labs(color = "Average Expression", size = "Percent Expressed") & ylab("Identity") & theme(axis.title.x = element_blank(), axis.text = element_text(size = 10)) & theme(legend.position = "right")
      p3 <- VlnPlot(scrna, features = i) & ggthemes::theme_hc() & scale_fill_manual(values = palette) & theme(legend.position = "right", axis.text = element_text(size = 10)) & labs(fill = "") & xlab("Identity") & ylab("Expression Level")

      # p1 <- plot_df %>% filter(gene == i) %>%       ggplot(aes(x=UMAP_1,y=UMAP_2,color=expr)) + geom_point(size=0.3) + theme_cellsnake_classic() + scale_color_continuous(type = "viridis") + labs(color="Expression") + theme(axis.text = element_text(size=12))
      # p2 <- plot_df %>% group_by(gene,clusters)  %>% summarise(percent=100*(length(expr[expr>0])/n()),average=mean(expr)) %>% filter(gene == i) %>%       ggplot(aes(x=gene,y=clusters,size=percent,color=average)) + geom_point() + theme_cellsnake_classic() +       scale_color_continuous(type = "viridis") + labs(color="Average Expression",size="Percent Expressed")+ ylab("Identity") + theme(axis.title.x = element_blank(),axis.text = element_text(size=12)) + theme(legend.position = "right")
      # p3 <- plot_df %>% ungroup() %>% filter(gene == i) %>%       ggplot(aes(x=clusters,y=log10(expr),fill=clusters)) + geom_violin() + theme_cellsnake_classic() + geom_jitter(size=0.3,shape=20, position=position_jitter(0.2)) +       scale_fill_manual(values = palette) + theme(legend.position = "right",axis.text = element_text(size = 12)) + labs(fill="") + xlab("Identity") + ylab("Log10-Expression Level")

      suppressWarnings(((p1 | p2) / p3) -> wp)
    },
    error = function(cond) {
      wp <- plot.new()
    },
    finally = {}
  ) -> wp
  ggsave(paste0(opt$output.plot.dir, "/", i, ".pdf"), wp, height = 5 + (n * 0.15), width = 7 + (n * 0.15), useDingbats = TRUE)
})