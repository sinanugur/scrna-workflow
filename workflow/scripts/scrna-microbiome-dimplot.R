#!/usr/bin/env Rscript


option_list <- list(
  optparse::make_option(c("--rds"),
    type = "character", default = NULL,
    help = "Processed rds file of a Seurat object", metavar = "character"
  ),
  optparse::make_option(c("--microbiome.rds"),
    type = "character", default = NULL,
    help = "Microbiome rds file", metavar = "character"
  ),
  optparse::make_option(c("--reduction.type"),
    type = "character", default = "umap",
    help = "Reduction type, umap or tsne", metavar = "character"
  ),
  optparse::make_option(c("--dimplot"),
    type = "character", default = "micdimplot.pdf",
    help = "Plot file name", metavar = "character"
  ),
  optparse::make_option(c("--tplot"),
    type = "character", default = "tplot.pdf",
    help = "Total microbiome dimplot file name", metavar = "character"
  ),
  optparse::make_option(c("--taxa"),
    type = "character", default = "genus",
    help = "Taxonomic level", metavar = "character"
  )
)



opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

if (is.null(opt$rds)) {
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file)", call. = FALSE)
}

require(patchwork)
require(Seurat)
# require(randomcoloR)

require(viridis)
require(tidyverse)

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

scrna <- readRDS(file = opt$rds)
DefaultAssay(scrna) <- "RNA"


microbiome <- readRDS(file = opt$microbiome.rds)



AddMetaData(scrna, microbiome %>% rownames_to_column("barcodes") %>% gather(taxa, umi, -barcodes) %>% dplyr::group_by(taxa) %>% dplyr::mutate(sum = sum(umi, na.rm = T)) %>% ungroup() %>%
  dplyr::mutate(taxa = ifelse(sum >= min(sort(unique(sum), decreasing = T)[1:11], na.rm = T), paste0(opt$taxa, "_", taxa), paste0(opt$taxa, "_", "others"))) %>%
  dplyr::select(-sum) %>% dplyr::group_by(barcodes, taxa) %>% dplyr::summarise(sum = sum(umi, na.rm = T)) %>% ungroup() %>% spread(taxa, sum) %>% column_to_rownames("barcodes")) -> scrna


# p1 <- DimPlot(scrna, reduction = opt$reduction.type, label = TRUE) & theme_cellsnake_classic() & scale_color_manual(values = palette)

scrna@meta.data %>%
  dplyr::mutate(`Total log2-Expression (Microbiome)` = log2(rowSums(across(starts_with(opt$taxa))) + 1)) %>%
  dplyr::mutate(across(contains("genus"), ~ replace(., .x == 0, NA))) -> scrna@meta.data


# scrna %>%
#  dplyr::select(barcodes = .cell, orig.ident, contains(opt$taxa), starts_with(opt$reduction.type)) %>%
#  gather(taxa, umi, starts_with(opt$taxa)) %>%
#  dplyr::select(barcodes, orig.ident, x = 3, y = 4, taxa, umi) %>%
#  replace(is.na(.), 0) %>%
#  ggplot(aes(x = x, y = y, color = log2(umi + 1))) +
#  geom_point(size = 0.2) +
#  labs(color = "Log2-UMI") +
#  theme(axis.text = element_text(size = 12)) +
#  scale_color_viridis(option = "magma", direction = -1, alpha = 0.8, na.value = "white") +
#  ggthemes::theme_few() +
#  facet_wrap(~taxa) -> p1


# ggsave(plot = p1, filename = opt$dimplot, width = 13, height = 9)




plotting_taxas <- scrna@meta.data %>%
  dplyr::select(starts_with(opt$taxa)) %>%
  select(
    where(
      ~ !all(is.na(.x))
    )
  ) %>%
  colnames() %>%
  unique()

print(plotting_taxas)

pdf(opt$dimplot, width = 7, height = 7)
for (i in plotting_taxas) {
  try({
    FeaturePlot(scrna, features = i, pt.size = 0.1, reduction = opt$reduction.type) &
      scale_color_continuous(type = "viridis", na.value = "gray96") -> p1

    p1 <- (p1 / guide_area()) + plot_layout(heights = c(2.5, 1), widths = c(1, 0.6), guides = "collect")
    print(p1)
  })
}
dev.off()



p2 <- FeaturePlot(scrna, features = "Total log2-Expression (Microbiome)", pt.size = 0.1, reduction = opt$reduction.type) &
  scale_color_continuous(type = "viridis", na.value = "gray96")

(p2 / guide_area()) + plot_layout(heights = c(2.5, 1), widths = c(1, 0.6), guides = "collect") -> p2

ggsave(plot = p2, filename = opt$tplot, width = 7.5, height = 8)