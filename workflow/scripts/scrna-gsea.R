#!/usr/bin/env Rscript

option_list <- list(
  optparse::make_option(c("--rds"),
    type = "character", default = NULL,
    help = "A list of RDS files of Seurat objects", metavar = "character"
  ),
  optparse::make_option(c("--gseafile"),
    type = "character", default = "c2.all.v2022.1.Hs.symbols.gmt",
    help = "GeneSetEnrichmentAnalysis file downloaded from http://www.gsea-msigdb.org/ ", metavar = "character"
  ),
  optparse::make_option(c("--idents"),
    type = "character", default = "seurat_clusters",
    help = "Groups for analysis", metavar = "character"
  ),
  optparse::make_option(c("--csv"),
    type = "character", default = NULL,
    help = "Celltypist prediction file", metavar = "character"
  ),
  optparse::make_option(c("--output.xlsx"),
    type = "character", default = NULL,
    help = "Output excel file name", metavar = "character"
  )
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

if (is.null(opt$rds)) {
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds and sampleid)", call. = FALSE)
}

require(cerebroApp)
require(openxlsx)
require(Seurat)
require(tidyverse)

scrna <- readRDS(file = opt$rds)

DefaultAssay(scrna) <- "RNA"

# Papers
# Diaz-Mejia, J. J. et al. Evaluation of methods to assign cell type labels to cell clusters from single-cell RNA-sequencing data. [version 3; peer review: 2 approved, 1 approved with reservations]. F1000Res. 8, ISCB Comm J-296 (2019).
# explanation: https://romanhaa.github.io/cerebroApp/reference/performGeneSetEnrichmentAnalysis.html

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


set.seed(155)
scrna <- performGeneSetEnrichmentAnalysis(
  object = scrna,
  GMT_file = opt$gseafile,
  groups = opt$idents,
  thresh_p_val = 0.05,
  name = "GSVA",
  thresh_q_val = 0.1
)

openxlsx::write.xlsx(scrna@misc[["enriched_pathways"]]$GSVA[opt$idents], opt$output.xlsx)