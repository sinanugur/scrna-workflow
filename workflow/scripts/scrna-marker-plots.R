#!/usr/bin/env Rscript
option_list = list(
  optparse::make_option(c("--resolution"), type="double", default=0.8, 
              help="Resolution [default= %default]", metavar="character"),

    optparse::make_option(c("--rds"), type="character", default=NULL, 
              help="Processed rds file of a Seurat object", metavar="character"),

    optparse::make_option(c("--sampleid"), type="character", default=NULL, 
              help="Sample ID", metavar="character"),
      optparse::make_option(c("--n"), type="integer", default=50, 
              help="How many features to plot per cluster [default= %default]", metavar="integer"),
    optparse::make_option(c("--xlsx"), type="character", default=NULL, 
              help="Excel table of markers", metavar="character")


)
 


opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$rds) || is.null(opt$sampleid) ){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file and sampleid)", call.=FALSE)
}

require(Seurat)
require(tidyverse)
require(viridis)
require(patchwork)




scrna=readRDS(file = opt$rds)

RNA_=paste0("RNA_snn_res.",opt$resolution)


Positive_Features=openxlsx::read.xlsx(opt$xlsx) %>% group_by(cluster) %>% slice_min(order_by = p_val_adj,n = opt$n) %>% select(cluster,gene) 


for(d in (Positive_Features %>% distinct(cluster) %>% pull())) {

    dir.create(paste0("results/",opt$sampleid,"/resolution-",opt$resolution,"/markers/","cluster",d,"/"),recursive=TRUE)
}

options(warn=-1)
suppressMessages(for (i in 1:nrow(Positive_Features)) {

    gene=Positive_Features[i,]$gene
    cluster=Positive_Features[i,]$cluster

p1 <- FeaturePlot(scrna, reduction = "umap", features=gene) + scale_color_continuous(type="viridis")
p2 <- DotPlot(scrna, features=gene)
p3 <- VlnPlot(scrna,features=gene)

suppressWarnings(((p1|p2)/p3) -> wp)

ggsave(paste0("results/",opt$sampleid,"/resolution-",opt$resolution,"/markers/cluster",cluster,"/",gene,".pdf"),wp,height=9,width=9)

})