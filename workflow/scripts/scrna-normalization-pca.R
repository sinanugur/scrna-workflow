#!/usr/bin/env Rscript


option_list <- list(
  optparse::make_option(c("--scale.factor"),
    type = "integer", default = 10000,
    help = "Scale factor [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--nfeatures"),
    type = "integer", default = 2000,
    help = "Highly variable features [default= %default]", metavar = "integer"
  ),
  optparse::make_option(c("--variable.selection.method"),
    type = "character", default = "vst",
    help = "Find variable features selection method [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--rds"),
    type = "character", default = NULL,
    help = "A RAW rds file of a Seurat object", metavar = "character"
  ),
  optparse::make_option(c("--normalization.method"),
    type = "character", default = "LogNormalize",
    help = "Normalization method[default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--doublet.filter"), action = "store_true", default = FALSE),
  optparse::make_option(c("--integration"), action = "store_true", default = FALSE),
  optparse::make_option(c("--umap"), action = "store_true", default = FALSE),
  optparse::make_option(c("--tsne"), action = "store_true", default = FALSE),
  optparse::make_option(c("--resolution"),
    type = "character", default = "0.8",
    help = "Resolution [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--output.rds"),
    type = "character", default = "output.rds",
    help = "Output RDS file name [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--cpu"),
    type = "integer", default = 5,
    help = "Number of CPU for parallel run [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--reference"),
    type = "character", default = "HumanPrimaryCellAtlasData",
    help = "SingleR reference", metavar = "character"
  )
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

if (is.null(opt$rds)) {
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file)", call. = FALSE)
}

require(optparse)
require(SingleR)
require(celldex)
require(tidyverse)
require(Seurat)
require(patchwork)



try(
  {
    source("workflow/scripts/scrna-functions.R")
  },
  silent = TRUE
)
try(
  {
    source(paste0(system("python -c 'import os; import cellsnake; print(os.path.dirname(cellsnake.__file__))'", intern = TRUE), "/scrna/workflow/scripts/scrna-functions.R"))
  },
  silent = TRUE
)


scrna <- readRDS(file = opt$rds)

if (isFALSE(opt$integration)) {
  scrna <- NormalizeData(scrna, normalization.method = opt$normalization.method, scale.factor = opt$scale.factor)
  scrna <- FindVariableFeatures(scrna, selection.method = opt$variable.selection.method, nfeatures = opt$nfeatures)
} else {
  try({
    DefaultAssay(scrna) <- "integrated"
  }) # for now only for Seurat, Harmony will come
}


# all.genes <- rownames(scrna) memory requirements can be large if using all genes
not.all.genes <- VariableFeatures(scrna) # only variable features

scrna <- ScaleData(scrna, features = not.all.genes)
scrna <- RunPCA(scrna, features = not.all.genes)
dimensionReduction <- function_pca_dimensions(scrna)
scrna <- FindNeighbors(scrna, dims = 1:dimensionReduction)

if (!opt$resolution %in% c("auto", "AUTO", "Auto")) {
  scrna <- FindClusters(scrna, resolution = as.numeric(opt$resolution))
} else {
  if (!requireNamespace("MultiKParallel", quietly = TRUE)) {
    remotes::install_github("sinanugur/MultiKParallel")
  }
  require(MultiKParallel)
  scrna_tmp <- FindClusters(scrna, resolution = seq(0.2, 2.5, 0.15))
  multik <- MultiKParallel(scrna_tmp, reps = 10, seed = 255, resolution = seq(0.2, 2.5, 0.15), numCores = opt$cpu, nPC = dimensionReduction)
  multik$k %>%
    tibble::as_tibble() %>%
    count(value) %>%
    filter(n == max(n)) %>%
    slice(1) %>%
    pull(value) -> K

  scrna_tmp[[]] %>%
    as.data.frame() %>%
    select(starts_with("RNA")) %>%
    rownames_to_column("barcode") %>%
    mutate(across(where(is.factor), as.character)) %>%
    as_tibble() %>%
    gather(res, clu, contains("RNA")) %>%
    distinct(res, clu) %>%
    count(res) %>%
    mutate(res = str_remove(res, "RNA_snn_res.")) %>%
    mutate(diff = abs(n - K)) %>%
    slice_min(order_by = diff) %>%
    pull(res) %>%
    as.numeric() -> optimal_resolution


  scrna <- FindClusters(scrna, resolution = optimal_resolution)
}


print("here")

if (opt$umap) {
  scrna <- RunUMAP(scrna, dims = 1:dimensionReduction)
}


if (opt$tsne) {
  scrna <- RunTSNE(scrna, dims = 1:dimensionReduction, check_duplicates = FALSE)
}



if (opt$doublet.filter) {
  if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
    remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
  }
  require(DoubletFinder)

  homotypic.prop <- modelHomotypic(scrna$seurat_clusters)
  nExp_poi <- round(0.075 * nrow(scrna@meta.data)) ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  scrna <- doubletFinder_v3(scrna, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  scrna <- doubletFinder_v3(scrna, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = paste0("pANN_0.25_0.09_", nExp_poi), sct = FALSE)


  scrna@meta.data %>%
    tibble::rownames_to_column("barcodes") %>%
    select(barcodes, starts_with("DF")) %>%
    select(barcodes, DoubletFinder = 2) -> doublet_df

  scrna@meta.data <- scrna@meta.data %>%
    select(!starts_with("DF")) %>%
    select(!starts_with("pANN")) %>%
    tibble::rownames_to_column("barcodes") %>%
    dplyr::left_join(doublet_df, by = "barcodes") %>%
    tibble::column_to_rownames("barcodes")


  subset(scrna, subset = DoubletFinder == "Singlet") -> scrna
  # scrna$DoubletFinder <- NULL
}




# RNA_=paste0("RNA_snn_res.",opt$resolution)

# metrics <- table(scrna@meta.data[["seurat_clusters"]], scrna@meta.data$orig.ident)

# p1 <- DimPlot(scrna, reduction = "pca", label = TRUE, label.size = 10)

# celltype annotation with SingleR
ref <- get(opt$reference)()



DefaultAssay(scrna) <- "RNA"
smObjSCE <- as.SingleCellExperiment(scrna)
pred <- SingleR(test = smObjSCE, ref = ref, labels = ref$label.fine)
AddMetaData(scrna, pred["pruned.labels"] %>% as.data.frame() %>% dplyr::select(singler = pruned.labels)) -> scrna

try({
  DefaultAssay(scrna) <- "integrated"
})


# output files
saveRDS(scrna, file = opt$output.rds)
# openxlsx::write.xlsx(metrics %>% as.data.frame() %>% select(Cluster = 1, everything()), file = opt$output.xlsx)
# ggsave(plot = p1, filename = opt$output.pca.plot, width = 9, height = 7)