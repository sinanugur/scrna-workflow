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
              help="Mitochondria filtering percentage [default= %default]", metavar="character")


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

 
#nFeature_RNA is the number of genes detected in each cell. nCount_RNA is the total number of molecules detected within a cell. 




try({scrna.data <- Read10X(data.dir = opt$data.dir)})
try({scrna.data <- Read10X_h5(data.dir = paste0(opt$data.dir,"filtered_feature_bc_matrix.h5"))})



scrna <- CreateSeuratObject(counts = scrna.data, project = opt$sampleid, min.cells = opt$min.cells, min.features = opt$min.features)


scrna[["percent.mt"]] <- PercentageFeatureSet(scrna, pattern = "^MT-")
scrna[["percent.rp"]] <- PercentageFeatureSet(scrna, pattern = "^RP[SL]")


output.dir=paste0("results/",opt$sampleid,"/technicals/")
dir.create(output.dir,recursive = T)

VlnPlot(scrna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4)


ggsave(paste0(output.dir,"before-qc-trimming-violinplot.pdf"), width = 10,height = 4)




 #minCov=opt$minCov #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
 #if(min(scrna$nCount_RNA)>=minCov){
 #   countLOW=min(scrna$nCount_RNA)
 # }else{
 #   countLOW=quantile(scrna$nCount_RNA, prob=0.01)  
 #}
 #countHIGH=quantile(scrna$nCount_RNA, prob=0.99) 
 #featureLOW=quantile(scrna$nFeature_RNA, prob=0.01)


lower_bound_nCount_RNA <- median(scrna$nCount_RNA) - 3 * mad(scrna$nCount_RNA, constant = 1)
upper_bound_nCount_RNA <- median(scrna$nCount_RNA) + 3 * mad(scrna$nCount_RNA, constant = 1)


lower_bound_nFeature_RNA <- median(scrna$nFeature_RNA) - 3 * mad(scrna$nFeature_RNA, constant = 1)
upper_bound_nFeature_RNA <- median(scrna$nFeature_RNA) + 3 * mad(scrna$nFeature_RNA, constant = 1)



##subset
#scrna <- subset(scrna, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & opt$percent.mt < opt$percent.mt)

scrna <- subset(scrna, subset = nFeature_RNA > lower_bound_nFeature_RNA & nFeature_RNA < upper_bound_nFeature_RNA & nCount_RNA > lower_bound_nCount_RNA  & nCount_RNA < upper_bound_nCount_RNA & percent.mt < opt$percent.mt)

VlnPlot(scrna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4)
ggsave(paste0(output.dir,"after-qc-trimming-violinplot.pdf"), width = 10,height = 4)


output.dir=paste0("analyses/raw/")
dir.create(output.dir,recursive = T)

saveRDS(scrna,file = paste0(output.dir,opt$sampleid,".rds"))
