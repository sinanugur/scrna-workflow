#!/usr/bin/env Rscript
option_list = list(

    optparse::make_option(c("--xlsx"), type="character", default=NULL, 
              help="Excel table of markers", metavar="character"),
    optparse::make_option(c("--output"), type="character", default=NULL, 
              help="Output excel file name", metavar="character")


)
 


opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$xlsx) || is.null(opt$output) ){
  optparse::print_help(opt_parser)
  stop("Arguments must be supplied", call.=FALSE)
}



require(Seurat)
require(tidyverse)
source("workflow/scripts/scrna-functions.R")





All_Features=openxlsx::read.xlsx(opt$xlsx)


All_Features %>% group_by(cluster) %>% split(.$cluster) %>% map(~function_gofunc(.)) %>% bind_rows(.id="cluster") -> go_enrichment

openxlsx::write.xlsx(go_enrichment,file=opt$output)






