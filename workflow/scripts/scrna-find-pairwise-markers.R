#!/usr/bin/env Rscript


option_list <- list(
    optparse::make_option(c("--rds"),
        type = "character", default = NULL,
        help = "Processed rds file of a Seurat object", metavar = "character"
    ),
    optparse::make_option(c("--logfc.threshold"),
        type = "double", default = 0.25,
        help = "LogFC [default= %default]", metavar = "character"
    ),
    optparse::make_option(c("--test.use"),
        type = "character", default = "wilcox",
        help = "Test use [default= %default]", metavar = "character"
    ),
    optparse::make_option(c("--output.rds"),
        type = "character", default = NULL,
        help = "RDS table of all markers", metavar = "character"
    ),
    optparse::make_option(c("--metadata.column"),
        type = "character", default = "condition",
        help = "Meta data column name for marker analysis", metavar = "character"
    ),
    optparse::make_option(c("--metadata"),
        type = "character", default = NULL,
        help = "Metadata filename", metavar = "character"
    )
)



opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

if (is.null(opt$rds)) {
    optparse::print_help(opt_parser)
    stop("At least one argument must be supplied (rds file)", call. = FALSE)
}

require(Seurat)
require(tidyverse)

function_read_metadata <- function(opt) {
    file_type <- tools::file_ext(opt$metadata)
    possible_separators <- c(",", ";", "|", "\t")

    if (tolower(file_type) %in% c("csv", "tsv", "txt")) {
        for (sep in possible_separators) {
            tryCatch(
                {
                    df <- read.csv(opt$metadata, sep = sep)
                    return(df)
                },
                error = function(e) {
                    message(paste("Unable to read file with separator", sep))
                }
            )
        }
    } else if (tolower(file_type) %in% c("xls", "xlsx")) {
        for (sep in possible_separators) {
            tryCatch(
                {
                    df <- openxlsx::read.xlsx(opt$metadata)
                    return(df)
                },
                error = function(e) {
                    message(paste("Unable to read file with separator", sep))
                }
            )
        }
    } else {
        message("Unsupported file type")
        return(NULL)
    }

    message("Unable to read file with any separator")
    return(NULL)
}


scrna <- readRDS(file = opt$rds)
DefaultAssay(scrna) <- "RNA"

function_read_metadata(opt) -> df


scrna@meta.data %>% dplyr::left_join(df, by = c("orig.ident" = names(df)[1])) -> scrna@meta.data


# RNA_=paste0("RNA_snn_res.",opt$resolution)
# Idents(object = scrna) <- scrna@meta.data[[RNA_]]
Idents(object = scrna) <- scrna@meta.data[[opt$metadata.column]]

combn(as.character(unique(Idents(scrna))), 2) -> pairwise



results <- list()
for (i in 1:ncol(pairwise)) {
    pair_ <- paste0(pairwise[1, i], " vs ", pairwise[2, i])
    markers <- FindMarkers(scrna, ident.1 = pairwise[1, i], ident.2 = pairwise[2, i], logfc.threshold = opt$logfc.threshold, test.use = opt$test.use)
    results[[pair_]] <- markers

    pair_ <- paste0(pairwise[2, i], " vs ", pairwise[1, i])
    markers <- FindMarkers(scrna, ident.1 = pairwise[2, i], ident.2 = pairwise[1, i], logfc.threshold = opt$logfc.threshold, test.use = opt$test.use)
    results[[pair_]] <- markers
}

results %>%
    map(~ rownames_to_column(., var = "gene")) %>%
    bind_rows(.id = "condition") -> all_markers


saveRDS(all_markers, file = opt$output.rds)