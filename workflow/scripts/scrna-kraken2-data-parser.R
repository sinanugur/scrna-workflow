#!/usr/bin/env Rscript

option_list <- list(
  optparse::make_option(c("--min.cells"),
    type = "integer", default = 1,
    help = "Min cells [default= %default]", metavar = "integer"
  ),
  optparse::make_option(c("--min.features"),
    type = "integer", default = 3,
    help = "Min features [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--data.dir"),
    type = "character", default = NULL,
    help = "Data directory", metavar = "character"
  ),
  optparse::make_option(c("--h5seurat"),
    type = "character", default = NULL,
    help = "Data directory", metavar = "character"
  ),
  optparse::make_option(c("--sampleid"),
    type = "character", default = NULL,
    help = "Sample ID", metavar = "character"
  ),
  optparse::make_option(c("--output.rds"),
    type = "character", default = "output.rds",
    help = "Output RDS file name [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--output.plot"),
    type = "character", default = "output.pdf",
    help = "Output barplot file name [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--taxa"),
    type = "character", default = "genus",
    help = "Taxonomic level", metavar = "character"
  )
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

if (is.null(opt$h5seurat)) {
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (data.dir and sampleid)", call. = FALSE)
}

require(optparse)
require(tidyverse)
require(Seurat)
require(SeuratDisk)
tryCatch(
  {
    source("workflow/scripts/scrna-functions.R")
  },
  error = function(cond) {
    source(paste0(system("python -c 'import os; import cellsnake; print(os.path.dirname(cellsnake.__file__))'", intern = TRUE), "/scrna/workflow/scripts/scrna-functions.R"))
  }
)
# nFeature_RNA is the number of genes detected in each cell. nCount_RNA is the total number of molecules detected within a cell.


# CreateSeuratObject(scrna.data,min.cells = 1,min.features = 5) -> scrna

CreateSeuratObject(LoadH5Seurat(opt$h5seurat)[["RNA"]]@counts, min.cells = opt$min.cells, min.features = opt$min.features) -> scrna

scrna <- RenameCells(object = scrna, add.cell.id = make.names(opt$sampleid)) # add cell.id to cell name


scrna[["RNA"]]@counts %>%
  as.matrix() %>%
  t() %>%
  as.data.frame() %>%
  select(-starts_with("Homo")) -> df

df %>%
  rownames_to_column("barcode") %>%
  gather(group, umi, -barcode) %>%
  group_by(group) %>%
  summarise(sum = log(sum(umi))) %>%
  arrange(desc(sum)) %>%
  slice_max(n = 50, order_by = sum) %>%
  ggplot(aes(reorder(group, sum), sum)) +
  geom_col() +
  coord_flip() +
  ggthemes::theme_few() +
  ylab("log-total UMI") +
  xlab(opt$taxa) +
  ggtitle(opt$sampleid) +
  theme(axis.title = element_text(size = 12)) -> p1

ggsave(plot = p1, filename = opt$output.plot, width = 6, height = 9)

saveRDS(df, file = opt$output.rds)