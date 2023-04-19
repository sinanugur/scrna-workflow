#!/usr/bin/env Rscript

option_list <- list(
  optparse::make_option(c("--min.cells"),
    type = "integer", default = 3,
    help = "Min cells [default= %default]", metavar = "integer"
  ),
  optparse::make_option(c("--min.features"),
    type = "integer", default = 200,
    help = "Min features, nFeature_RNA [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--max.features"),
    type = "integer", default = Inf,
    help = "Max features, nFeature_RNA [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--max.molecules"),
    type = "integer", default = Inf,
    help = "Max molecules, nCount_RNA [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--min.molecules"),
    type = "integer", default = 0,
    help = "Min molecules, nCount_RNA [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--data.dir"),
    type = "character", default = NULL,
    help = "Data directory", metavar = "character"
  ),
  optparse::make_option(c("--sampleid"),
    type = "character", default = NULL,
    help = "Sample ID", metavar = "character"
  ),
  optparse::make_option(c("--percent.mt"),
    type = "character", default = "10",
    help = "Max mitochondrial gene percentage, smaller than [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--percent.rp"),
    type = "double", default = 0,
    help = "Min ribosomal gene percentage, greater than [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--before.violin.plot"),
    type = "character", default = "before.violin.pdf",
    help = "Violin plot name [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--after.violin.plot"),
    type = "character", default = "after.violin.pdf",
    help = "Violin plot name [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--output.rds"),
    type = "character", default = "output.rds",
    help = "Output RDS file name [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--plot.mtplot"),
    type = "character", default = "plot.mtplot.pdf",
    help = "Violin plot name [default= %default]", metavar = "character"
  )
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

if (is.null(opt$data.dir) || is.null(opt$sampleid)) {
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (data.dir and sampleid)", call. = FALSE)
}

require(optparse)
require(tidyverse)
require(Seurat)
require(patchwork)
require(tools)
require(data.table)


# nFeature_RNA is the number of genes detected in each cell. nCount_RNA is the total number of molecules detected within a cell.



function_read_input <- function(opt) {
  try(
    {
      scrna.data <- Read10X(data.dir = opt$data.dir)
      return(scrna.data)
    },
    silent = TRUE
  )
  try(
    {
      scrna.data <- Read10X_h5(filename = paste0(opt$data.dir, "/filtered_feature_bc_matrix.h5"))
      return(scrna.data)
    },
    silent = TRUE
  )
  try(
    {
      scrna.data <- Read10X_h5(filename = opt$data.dir)
      return(scrna.data)
    },
    silent = TRUE
  )

  try(
    {
      x <- tolower(file_ext(opt$data.dir))

      if (x %in% c("h5")) {
        scrna.data <- Read10X_h5(filename = opt$data.dir)
      } else if (x %in% c("gz")) {
        scrna.data <- fread(cmd = paste("gunzip -dc", opt$data.dir))

        scrna.data <- scrna.data %>% column_to_rownames("V1")
      } else if (x %in% c("zip")) {
        scrna.data <- fread(cmd = paste("unzip -p", opt$data.dir))

        scrna.data <- scrna.data %>% column_to_rownames("V1")
      } else if (x %in% c("csv", "tsv")) {
        scrna.data <- fread(paste(opt$data.dir))

        scrna.data <- scrna.data %>% column_to_rownames("V1")
      }
      return(scrna.data)
    },
    silent = TRUE
  )
}



function_read_input(opt) -> scrna.data




scrna <- CreateSeuratObject(counts = scrna.data, project = make.names(opt$sampleid), min.cells = opt$min.cells, min.features = opt$min.features)
rm(scrna.data)


scrna <- RenameCells(object = scrna, add.cell.id = make.names(opt$sampleid)) # add cell.id to cell name


scrna[["percent.mt"]] <- PercentageFeatureSet(scrna, pattern = "^[Mm][Tt]-")
scrna[["percent.rp"]] <- PercentageFeatureSet(scrna, pattern = "(?i)(^RP[SL])")

VlnPlot(scrna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol = 4)


ggsave(opt$before.violin.plot, width = 10, height = 4)

# auto detection is obsolote but keep it for later use
# lower_bound_nCount_RNA <- median(scrna$nCount_RNA) - 3 * mad(scrna$nCount_RNA, constant = 1)
# upper_bound_nCount_RNA <- median(scrna$nCount_RNA) + 3 * mad(scrna$nCount_RNA, constant = 1)
# lower_bound_nFeature_RNA <- median(scrna$nFeature_RNA) - 3 * mad(scrna$nFeature_RNA, constant = 1)
# upper_bound_nFeature_RNA <- median(scrna$nFeature_RNA) + 3 * mad(scrna$nFeature_RNA, constant = 1)



## subset
scrna <- subset(scrna, subset = nFeature_RNA < opt$max.features & nFeature_RNA >= opt$min.features & nCount_RNA < opt$max.molecules & nCount_RNA >= opt$min.molecules)

# scrna <- subset(scrna, subset = nFeature_RNA > lower_bound_nFeature_RNA & nFeature_RNA < upper_bound_nFeature_RNA & nCount_RNA > lower_bound_nCount_RNA  & nCount_RNA < upper_bound_nCount_RNA & percent.mt < opt$percent.mt)

if (opt$percent.mt %in% c("auto", "Auto", "AUTO")) {
  if (isFALSE(all(scrna@meta.data$percent.mt == 0))) {
    require(SingleCellExperiment)
    require(miQC)
    require(scater)


    smObjSCE <- as.SingleCellExperiment(scrna)
    mt_genes <- grepl("^[Mm][Tt]-", rownames(smObjSCE))
    feature_ctrls <- list(mito = rownames(smObjSCE)[mt_genes])
    smObjSCE <- addPerCellQC(smObjSCE, subsets = feature_ctrls)

    tryCatch(
      {
        model <- mixtureModel(smObjSCE)
        p1 <- plotModel(smObjSCE, model)
        p2 <- plotMetrics(smObjSCE)
        ggsave(filename = opt$plot.mtplot, p1 + p2, width = 10, height = 4)
        smObjSCE <- filterCells(smObjSCE, model)
        scrna <- scrna[, colnames(smObjSCE)]
        return(scrna)
      },
      error = function(a) {
        upper_bound_MT <- median(scrna$percent.mt) + 1 * mad(scrna$percent.mt, constant = 1) # miQC failed, use median absolute deviation

        scrna <- subset(scrna, subset = percent.mt <= upper_bound_MT)
        p1 <- plot.new()
        p2 <- plotMetrics(smObjSCE)

        ggsave(filename = opt$plot.mtplot, p1 + p2, width = 10, height = 4)
        return(scrna)
      }
    ) -> scrna
  }
} else {
  scrna <- subset(scrna, subset = percent.mt <= as.numeric(opt$percent.mt))
}

scrna <- subset(scrna, subset = percent.rp >= opt$percent.rp)



VlnPlot(scrna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol = 4)

ggsave(opt$after.violin.plot, width = 10, height = 4)


saveRDS(scrna, file = opt$output.rds)