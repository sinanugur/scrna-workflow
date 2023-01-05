#!/usr/bin/env Rscript

option_list = list(
    optparse::make_option(c("--rds"), type="character", default=NULL, 
              help="A list of RDS files of Seurat objects", metavar="character"),
    optparse::make_option(c("--gseafile"), type="character", default="c2.all.v2022.1.Hs.symbols.gmt", 
              help="GeneSetEnrichmentAnalysis file downloaded from http://www.gsea-msigdb.org/ ", metavar="character"),
    optparse::make_option(c("--group"), type="character", default="seurat_clusters", 
              help="Groups for analysis", metavar="character"),
        optparse::make_option(c("--output.dir"), type="character", default=NULL, 
              help="Output excel file name", metavar="character")
)
 
opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$rds)) {
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds and sampleid)", call.=FALSE)
}

require(cerebroApp)
require(openxlsx)
require(Seurat)
require(tidyverse)
try({source("workflow/scripts/scrna-functions.R")},silent=TRUE)
try({source(paste0(system("python -c 'import os; import cellsnake; print(os.path.dirname(cellsnake.__file__))'", intern = TRUE),"/scrna/workflow/scripts/scrna-functions.R"))})


scrna=readRDS(file = opt$rds)

#Papers
#Diaz-Mejia, J. J. et al. Evaluation of methods to assign cell type labels to cell clusters from single-cell RNA-sequencing data. [version 3; peer review: 2 approved, 1 approved with reservations]. F1000Res. 8, ISCB Comm J-296 (2019).
#explanation: https://romanhaa.github.io/cerebroApp/reference/performGeneSetEnrichmentAnalysis.html


set.seed(155)
scrna <- performGeneSetEnrichmentAnalysis(
  object = scrna,
  GMT_file = opt$gseafile,
  groups = opt$group,
  thresh_p_val = 0.05,
  name = "GSVA",
  thresh_q_val = 0.1
)

dir.create(opt$output.dir,recursive = TRUE)
write.xlsx(scrna@misc[["enriched_pathways"]]$GSVA[opt$group], paste0(opt$output.dir ,"/gsea-", opt$group, "-output.xlsx"))

