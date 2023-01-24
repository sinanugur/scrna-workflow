#!/usr/bin/env Rscript


option_list = list(
  optparse::make_option(c("--resolution"), type="double", default=0.8, 
              help="Resolution [default= %default]", metavar="character"),
    optparse::make_option(c("--rds"), type="character", default=NULL, 
              help="Processed rds file of a Seurat object", metavar="character"),
        optparse::make_option(c("--reduction.type"), type="character", default="umap", 
              help="Reduction type, umap or tsne", metavar="character"),
      optparse::make_option(c("--output.reduction.plot"), type="character", default="reduction.pdf", 
              help="Plot file name", metavar="character")


)
 


opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$rds)){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file)", call.=FALSE)
}

require(Seurat)
require(tidyverse)
require(randomcoloR)

try({source("workflow/scripts/scrna-functions.R")},silent=TRUE)
try({source(paste0(system("python -c 'import os; import cellsnake; print(os.path.dirname(cellsnake.__file__))'", intern = TRUE),"/scrna/workflow/scripts/scrna-functions.R"))},silent=TRUE)

scrna=readRDS(file = opt$rds)
DefaultAssay(scrna) <- "RNA"

n<-length(Idents(scrna) %>% unique())
set.seed(149)
palette <- sort(distinctColorPalette(n))


names(palette)=Idents(scrna) %>% unique() %>% sort() %>% as.character()
print(palette)

p1 <- DimPlot(scrna, reduction = opt$reduction.type, label = TRUE) & theme_cellsnake_classic() & scale_color_manual(values = palette) 




ggsave(plot =p1,filename=opt$output.reduction.plot,width=9,height=7)

