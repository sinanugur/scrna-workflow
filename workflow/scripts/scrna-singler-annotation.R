#!/usr/bin/env Rscript

option_list <- list(
    optparse::make_option(c("--rds"),
        type = "character", default = NULL,
        help = "A list of RDS files of Seurat objects", metavar = "character"
    ),
    optparse::make_option(c("--output"),
        type = "character", default = "pred.rds",
        help = "Output prediction file", metavar = "character"
    ),
    optparse::make_option(c("--reference"),
        type = "character", default = "HumanPrimaryCellAtlasData",
        help = "SingleR reference", metavar = "character"
    ),
    optparse::make_option(c("--granulation"),
        type = "character", default = "label.main",
        help = "SingleR granulation level", metavar = "character"
    )
)


opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

if (is.null(opt$rds)) {
    optparse::print_help(opt_parser)
    stop("At least one argument must be supplied (rds and sampleid)", call. = FALSE)
}


require(optparse)
require(SingleR)
require(SingleCellExperiment)
require(Seurat)
scrna <- readRDS(file = opt$rds)
DefaultAssay(scrna) <- "RNA"



# celltype annotation with SingleR
ref <- get(opt$reference)()

smObjSCE <- as.SingleCellExperiment(scrna)
pred <- SingleR(test = smObjSCE, ref = ref, labels = ref[[opt$granulation]])

saveRDS(pred, file = opt$output)