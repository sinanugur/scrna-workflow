#!/usr/bin/env Rscript


option_list = list(
  optparse::make_option(c("--min.cells"), type="integer", default=3, 
              help="Min cells [default= %default]", metavar="integer"),
    optparse::make_option(c("--min.features"), type="integer", default=200, 
              help="Min features [default= %default]", metavar="character"),
    optparse::make_option(c("--data.dir"), type="character", default=NULL, 
              help="Data directory", metavar="character"),
  optparse::make_option(c("--scale.factor"), type="integer", default=10000, 
              help="Scale factor [default= %default]", metavar="character"),
    optparse::make_option(c("--cpu"), type="integer", default=5, 
              help="Number of CPU for parallel run [default= %default]", metavar="character"),
  optparse::make_option(c("--nfeatures"), type="integer", default=2000, 
              help="Highly variable features [default= %default]", metavar="integer"),
    optparse::make_option(c("--rds"), type="character", default=NULL, 
              help="A RAW rds file of a Seurat object", metavar="character"),
    optparse::make_option(c("--normalization.method"), type="character", default="LogNormalize", 
              help="Normalization method[default= %default]", metavar="character"),
    optparse::make_option(c("--doublet.filter"), action = "store_true", default = FALSE),
    optparse::make_option(c("--sampleid"), type="character", default=NULL, 
              help="Sample ID", metavar="character"),
    optparse::make_option(c("--output.tsv"), type="character", default=NULL, 
              help="Output TSV file name", metavar="character"),
    optparse::make_option(c("--percent.mt"), type="double", default=10, 
              help="Mitochondria filtering percentage [default= %default]", metavar="character")



)
 
opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$data.dir) ){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file)", call.=FALSE)
}

require(optparse)
require(Seurat)
require(patchwork)
require(DoubletFinder)
require(MultiKParallel)
require(tidyverse)
try({source("workflow/scripts/scrna-functions.R")})
try({source(paste0(system("python -c 'import os; import cellsnake; print(os.path.dirname(cellsnake.__file__))'", intern = TRUE),"/scrna/workflow/scripts/scrna-functions.R"))})


try({scrna.data <- Read10X(data.dir = opt$data.dir)})
try({scrna.data <- Read10X_h5(filename = paste0(opt$data.dir,"/filtered_feature_bc_matrix.h5"))})



scrna <- CreateSeuratObject(counts = scrna.data, project = opt$sampleid, min.cells = opt$min.cells, min.features = opt$min.features)
scrna <- NormalizeData(scrna, normalization.method = opt$normalization.method, scale.factor = opt$scale.factor)
scrna <- FindVariableFeatures(scrna, selection.method = "vst", nfeatures = opt$nfeatures)

all.genes <- rownames(scrna)
scrna <- ScaleData(scrna, features = all.genes)
scrna <- RunPCA(scrna, features = VariableFeatures(object = scrna))
dimensionReduction=function_pca_dimensions(scrna)
scrna <- FindNeighbors(scrna, dims = 1:dimensionReduction)
scrna <- FindClusters(scrna, resolution = seq(0.2,2.5,0.15))


multik <- MultiKParallel(scrna, reps=10,seed = 255,resolution = seq(0.2,2.5,0.15),numCores = opt$cpu,nPC = dimensionReduction)


multik$k %>% tibble::as_tibble() %>% count(value) %>% filter(n == max(n)) %>% slice(1) %>% pull(value) -> K




scrna[[]] %>% as.data.frame() %>% select(starts_with("RNA")) %>% rownames_to_column("barcode") %>% mutate(across(where(is.factor),as.character)) %>% as_tibble() %>%
 gather(res,clu,contains("RNA")) %>% distinct(res,clu) %>% count(res) %>% mutate(res=str_remove(res,"RNA_snn_res.")) %>% mutate(diff= abs(n - K)) %>% slice_min(order_by = diff) %>% pull(res) %>% as.numeric() -> optimal_resolution



write_tsv(x = data.frame(sample_name=opt$sampleid, resolution=optimal_resolution,MT=opt$percent.mt),file = opt$output.tsv)