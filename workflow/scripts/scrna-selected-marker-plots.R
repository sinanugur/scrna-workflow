#!/usr/bin/env Rscript
option_list = list(
  optparse::make_option(c("--resolution"), type="double", default=0.8, 
              help="Resolution [default= %default]", metavar="character"),

    optparse::make_option(c("--rds"), type="character", default=NULL, 
              help="Processed rds file of a Seurat object", metavar="character"),

    optparse::make_option(c("--sampleid"), type="character", default=NULL, 
              help="Sample ID", metavar="character")


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

domanska_muscularis=data.frame(gene=c("FCER1A","CDC1C","CLEC10A","CCL3L1","CCL3","CCL4L2","MT1X","MT1E","CTSL","RGS1","FOS","APOE","DNASE1L3","MMP9","LYZ","AREG","EREG","CCL20","S100A9","S100A8","EREG","FCN1","VCAN","LYZ","HSPA1A","HSPA6","HSPA1B","LYVE1","MARCO","COLEC12"),group=c(3,3,3,5,5,5,6,6,6,7,7,7,4,4,4,1,1,1,0,0,0,2,2,2,9,9,9,11,11,11))
domanska_mucosas=data.frame(gene=c("LYVE1","F13A1","FOLR2","SELENOP","APOE","SLC40A1","C1QB","DAB2","PDK4","SPP1","ACP5","CD9","FCER1A","CD1C","CLEC10A","HSPA6","DNAJB1","HSPA1B","S100A8","S100A9","S100A12","EREG","G0S2","FCN1","CCL20","IL1B","IL23A","CXCL10","CXCL9","GBP1"),group=c(12,12,12,11,11,11,10,10,10,9,9,9,8,8,8,5,5,5,0,0,0,3,3,3,2,2,2,6,6,6))

domanska_markers=bind_rows(domanska_mucosas,domanska_muscularis) %>% distinct(gene) %>% pull()


scrna=readRDS(file = opt$rds)

RNA_=paste0("RNA_snn_res.",opt$resolution)

output.dir=paste0("results/",opt$sampleid,"/resolution-",opt$resolution,"/")
dir.create(paste0(output.dir,"selected-markers/plots/"),recursive = T)

Idents(object = scrna) <- scrna@meta.data[[RNA_]]



options(warn=-1)
suppressMessages(for (i in domanska_markers) {


try({

p1 <- FeaturePlot(scrna, reduction = "umap", features=i) + scale_color_continuous(type="viridis")
p2 <- DotPlot(scrna, features=i)
p3 <- VlnPlot(scrna,features=i)

suppressWarnings(((p1|p2)/p3) -> wp)

ggsave(paste0("results/",opt$sampleid,"/resolution-",opt$resolution,"/selected-markers/plots/",i,".pdf"),wp,height=9,width=9)

})

})

domanska_markers=sort(intersect(domanska_markers,rownames(scrna)))

DotPlot(scrna, features = domanska_markers, dot.scale = 8) + RotatedAxis()
ggsave(paste0(output.dir,"/selected-markers/","selected-markers-dotplot.pdf"),width = 13,height = 8)
