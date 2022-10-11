#!/usr/bin/env Rscript
option_list = list(
  optparse::make_option(c("--resolution"), type="double", default=0.8, 
              help="Resolution [default= %default]", metavar="character"),
    optparse::make_option(c("--rds"), type="character", default=NULL, 
              help="Processed rds file of a Seurat object", metavar="character"),
    
    optparse::make_option(c("--tsv"), type="character", default=NULL, 
              help="A text file contains the gene list", metavar="character"),
    
    optparse::make_option(c("--output.dotplot"), type="character", default=NULL, 
              help="Output dotplot name", metavar="character"),
          
    optparse::make_option(c("--output.plot.dir"), type="character", default=NULL, 
              help="Output plot directory", metavar="character"),
            optparse::make_option(c("--reduction.type"), type="character", default="umap", 
              help="Reduction type, umap or tsne", metavar="character")


)
 


opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$rds)){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file)", call.=FALSE)
}

require(Seurat)
require(tidyverse)
require(viridis)

markers=read_tsv(opt$tsv,col_names=FALSE) %>% pull()


scrna=readRDS(file = opt$rds)

RNA_=paste0("RNA_snn_res.",opt$resolution)


Idents(object = scrna) <- scrna@meta.data[[RNA_]]



options(warn=-1)
suppressMessages(for (i in markers) {

dir.create(opt$output.plot.dir,recursive = TRUE)

try({

p1 <- FeaturePlot(scrna, reduction = opt$reduction.type, features=i) + scale_color_continuous(type="viridis")
p2 <- DotPlot(scrna, features=i)
p3 <- VlnPlot(scrna,features=i)

suppressWarnings(((p1|p2)/p3) -> wp)

ggsave(paste0(opt$output.plot.dir,"/",i,".",opt$reduction.type,".pdf"),wp,height=9,width=9)

})

})

markers=sort(intersect(markers,rownames(scrna)))

DotPlot(scrna, features = markers, dot.scale = 8) + RotatedAxis()
ggsave(opt$output.dotplot,width = 13,height = 8)
