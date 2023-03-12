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
  optparse::make_option(c("--sigplot"),
    type = "character", default = "sigplot.pdf",
    help = "Plot file name", metavar = "character"
  ),
  optparse::make_option(c("--sigtable"),
    type = "character", default = NULL,
    help = "Excel table of positive markers", metavar = "character"
  ),
  optparse::make_option(c("--taxa"),
    type = "character", default = "genus",
    help = "Taxonomic level", metavar = "character"
  ),
  optparse::make_option(c("--idents"),
    type = "character", default = "seurat_clusters",
    help = "Meta data column name", metavar = "character"
  )
)



opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

if (is.null(opt$rds)) {
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file)", call. = FALSE)
}

require(Seurat)
# require(randomcoloR)
require(tidyseurat)
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

scrna %>%
  dplyr::select(one_of(opt$idents), starts_with(opt$taxa)) %>%
  gather(taxa, umi, starts_with(opt$taxa)) %>%
  group_by(across(opt$idents), taxa) %>%
  dplyr::mutate(total = sum(umi, na.rm = T)) %>%
  group_by(taxa, across(opt$idents)) %>%
  dplyr::mutate(cell = n()) %>%
  dplyr::ungroup() %>%
  distinct(across(opt$idents), taxa, total, cell) %>%
  group_by(taxa) %>%
  dplyr::mutate(v3 = sum(total) - total, v4 = sum(cell) - cell) %>%
  rowwise() %>%
  dplyr::mutate(p = fisher.test(matrix(c(total, cell, v3, v4), ncol = 2), alternative = "greater")$p.value) %>%
  ungroup() %>%
  dplyr::mutate(p = p.adjust(p)) %>%
  arrange(p) -> df


df %>% ggplot(aes(x = get(opt$idents), y = -log10(p + 1e-200))) +
  geom_col() +
  geom_hline(yintercept = -log10(0.05), color = "red") +
  facet_wrap(~taxa) +
  theme_bw() +
  coord_flip() +
  ylab("Adjusted log P value") +
  xlab("Identity") -> p1

ggsave(plot = p1, filename = opt$sigplot, width = 12, height = 8)


openxlsx::write.xlsx(df, file = opt$sigtable)