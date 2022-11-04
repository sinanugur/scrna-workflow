#!/usr/bin/env Rscript

option_list = list(
  optparse::make_option(c("--min.cells"), type="integer", default=3, 
              help="Min cells [default= %default]", metavar="integer"),
    optparse::make_option(c("--min.features"), type="integer", default=200, 
              help="Min features [default= %default]", metavar="character"),
    optparse::make_option(c("--data.dir"), type="character", default=NULL, 
              help="Data directory", metavar="character"),
    optparse::make_option(c("--sampleid"), type="character", default=NULL, 
              help="Sample ID", metavar="character"),
    optparse::make_option(c("--percent.mt"), type="double", default=10, 
              help="Mitochondria filtering percentage, smaller than [default= %default]", metavar="character"),
    optparse::make_option(c("--percent.rp"), type="double", default=0, 
              help="Ribosomal filtering percentage, greater than [default= %default]", metavar="character"),

    optparse::make_option(c("--before.violin.plot"), type="character", default="before.violin.pdf", 
              help="Violin plot name [default= %default]", metavar="character"),
    optparse::make_option(c("--after.violin.plot"), type="character", default="after.violin.pdf", 
              help="Violin plot name [default= %default]", metavar="character")
              ,
     optparse::make_option(c("--output.rds"), type="character", default="output.rds", 
              help="Output RDS file name [default= %default]", metavar="character"),

    optparse::make_option(c("--auto.mt.filter"), action="store_true", default=FALSE,
              help="Should the program do auto mitochondria filtering [default %default]"),

    optparse::make_option(c("--plot.mtplot"), type="character", default="plot.mtplot.pdf", 
              help="Violin plot name [default= %default]", metavar="character")


)
 
opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$data.dir) || is.null(opt$sampleid) ){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (data.dir and sampleid)", call.=FALSE)
}

require(optparse)
require(tidyverse)
require(Seurat)
require(patchwork)
require(SingleCellExperiment)
require(miQC)
require(scater)

 
#nFeature_RNA is the number of genes detected in each cell. nCount_RNA is the total number of molecules detected within a cell. 




try({scrna.data <- Read10X(data.dir = opt$data.dir)})
try({scrna.data <- Read10X_h5(filename = paste0(opt$data.dir,"/filtered_feature_bc_matrix.h5"))})



scrna <- CreateSeuratObject(counts = scrna.data, project = opt$sampleid, min.cells = opt$min.cells, min.features = opt$min.features)


scrna[["percent.mt"]] <- PercentageFeatureSet(scrna, pattern = "^MT-")
scrna[["percent.rp"]] <- PercentageFeatureSet(scrna, pattern = "^RP[SL]")

VlnPlot(scrna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4)


ggsave(opt$before.violin.plot, width = 10,height = 4)

#auto detection is obsolote but keep it for later use
#lower_bound_nCount_RNA <- median(scrna$nCount_RNA) - 3 * mad(scrna$nCount_RNA, constant = 1)
#upper_bound_nCount_RNA <- median(scrna$nCount_RNA) + 3 * mad(scrna$nCount_RNA, constant = 1)
#lower_bound_nFeature_RNA <- median(scrna$nFeature_RNA) - 3 * mad(scrna$nFeature_RNA, constant = 1)
#upper_bound_nFeature_RNA <- median(scrna$nFeature_RNA) + 3 * mad(scrna$nFeature_RNA, constant = 1)



##subset
#scrna <- subset(scrna, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & opt$percent.mt < opt$percent.mt)

#scrna <- subset(scrna, subset = nFeature_RNA > lower_bound_nFeature_RNA & nFeature_RNA < upper_bound_nFeature_RNA & nCount_RNA > lower_bound_nCount_RNA  & nCount_RNA < upper_bound_nCount_RNA & percent.mt < opt$percent.mt)

if (opt$auto.mt.filter) {

 smObjSCE = as.SingleCellExperiment(scrna)
  mt_genes <- grepl("^MT-",  rownames(smObjSCE))
  feature_ctrls <- list(mito = rownames(smObjSCE)[mt_genes])
  smObjSCE <- addPerCellQC(smObjSCE, subsets = feature_ctrls)


model <- mixtureModel(smObjSCE)

p1 <- plotModel(smObjSCE, model) 
p2 <- plotMetrics(smObjSCE)
ggsave(filename=opt$plot.mtplot, p1 + p2,width = 10,height = 4)


smObjSCE <- filterCells(smObjSCE, model)

scrna <- scrna[,colnames(smObjSCE)]


scrna <- subset(scrna, subset = percent.rp > opt$percent.rp)

} else {


scrna <- subset(scrna, subset = percent.mt < opt$percent.mt & percent.rp > opt$percent.rp)


}



VlnPlot(scrna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4)

ggsave(opt$after.violin.plot, width = 10,height = 4)


saveRDS(scrna,file = opt$output.rds)







