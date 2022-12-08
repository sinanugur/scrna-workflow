#!/usr/bin/env Rscript


option_list = list(
  optparse::make_option(c("--scale.factor"), type="integer", default=10000, 
              help="Scale factor [default= %default]", metavar="character"),
  optparse::make_option(c("--nfeatures"), type="integer", default=2000, 
              help="Highly variable features [default= %default]", metavar="integer"),
    optparse::make_option(c("--rds"), type="character", default=NULL, 
              help="A RAW rds file of a Seurat object", metavar="character"),
    optparse::make_option(c("--normalization.method"), type="character", default="LogNormalize", 
              help="Normalization method[default= %default]", metavar="character"),
    optparse::make_option(c("--doublet.filter"), action = "store_true", default = FALSE),
    optparse::make_option(c("--umap"), action = "store_true", default = FALSE),
    optparse::make_option(c("--tsne"), action = "store_true", default = FALSE),
    optparse::make_option(c("--resolution"), type="character", default="0.8", 
              help="Resolution [default= %default]", metavar="character"),
    optparse::make_option(c("--output.rds"), type="character", default="output.rds", 
              help="Output RDS file name [default= %default]", metavar="character"),
    optparse::make_option(c("--output.xlsx"), type="character", default=NULL, 
              help="Excel table of markers", metavar="character"),
    optparse::make_option(c("--output.pca.plot"), type="character", default="pca.pdf", 
              help="PCA plot file name", metavar="character"),
    optparse::make_option(c("--cpu"), type="integer", default=5, 
              help="Number of CPU for parallel run [default= %default]", metavar="character")


)
 
opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$rds) ){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file)", call.=FALSE)
}

require(optparse)
require(tidyverse)
require(Seurat)
require(patchwork)
require(DoubletFinder)
try({source("workflow/scripts/scrna-functions.R")})
try({source(paste0(system("python -c 'import os; import cellsnake; print(os.path.dirname(cellsnake.__file__))'", intern = TRUE),"/scrna/workflow/scripts/scrna-functions.R"))})


scrna=readRDS(file = opt$rds)


scrna <- NormalizeData(scrna, normalization.method = opt$normalization.method, scale.factor = opt$scale.factor)
scrna <- FindVariableFeatures(scrna, selection.method = "vst", nfeatures = opt$nfeatures)

all.genes <- rownames(scrna)
scrna <- ScaleData(scrna, features = all.genes)
scrna <- RunPCA(scrna, features = VariableFeatures(object = scrna))
dimensionReduction=function_pca_dimensions(scrna)
scrna <- FindNeighbors(scrna, dims = 1:dimensionReduction)

if(opt$resolution != "auto") {
scrna <- FindClusters(scrna, resolution = as.numeric(opt$resolution))
} else {
require(MultiKParallel)
scrna_tmp <- FindClusters(scrna, resolution = seq(0.2,2.5,0.15))
multik <- MultiKParallel(scrna_tmp, reps=10,seed = 255,resolution = seq(0.2,2.5,0.15),numCores = opt$cpu,nPC = dimensionReduction)
multik$k %>% tibble::as_tibble() %>% count(value) %>% filter(n == max(n)) %>% slice(1) %>% pull(value) -> K

scrna_tmp[[]] %>% as.data.frame() %>% select(starts_with("RNA")) %>% rownames_to_column("barcode") %>% mutate(across(where(is.factor),as.character)) %>% as_tibble() %>%
 gather(res,clu,contains("RNA")) %>% distinct(res,clu) %>% count(res) %>% mutate(res=str_remove(res,"RNA_snn_res.")) %>% mutate(diff= abs(n - K)) %>% slice_min(order_by = diff) %>% pull(res) %>% as.numeric() -> optimal_resolution


scrna <- FindClusters(scrna, resolution = optimal_resolution)

}


if(opt$umap) {scrna <- RunUMAP(scrna, dims = 1:dimensionReduction)}
if(opt$tsne) {scrna <- RunTSNE(scrna, dims = 1:dimensionReduction)}




homotypic.prop <- modelHomotypic(scrna$seurat_clusters)         
nExp_poi <- round(0.075*nrow(scrna@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
scrna <- doubletFinder_v3(scrna, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
scrna <- doubletFinder_v3(scrna, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = paste0("pANN_0.25_0.09_",nExp_poi), sct = FALSE)


scrna@meta.data %>% tibble::rownames_to_column("barcodes") %>% select(barcodes,starts_with("DF")) %>% select(barcodes,DoubletFinder=2) -> Doublet_Df

scrna@meta.data <- scrna@meta.data %>% select(!starts_with("DF")) %>% select(!starts_with("pANN")) %>%
  tibble::rownames_to_column("barcodes") %>%
  dplyr::left_join(Doublet_Df, by = "barcodes") %>%
  tibble::column_to_rownames("barcodes")

if(opt$doublet.filter) {
subset(scrna, subset=DoubletFinder == "Singlet") -> scrna
#scrna$DoubletFinder <- NULL

}




#RNA_=paste0("RNA_snn_res.",opt$resolution)

metrics=table(scrna@meta.data[["seurat_clusters"]], scrna@meta.data$orig.ident)

p1 <- DimPlot(scrna, reduction = "pca", label = TRUE, label.size = 10) 


#output files
saveRDS(scrna,file = opt$output.rds)
openxlsx::write.xlsx(metrics %>% as.data.frame() %>% select(Cluster=1,everything()),file=opt$output.xlsx)
ggsave(plot =p1,filename=opt$output.pca.plot,width=9,height=7)



