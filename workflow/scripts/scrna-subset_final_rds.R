#!/usr/bin/env Rscript


option_list <- list(
    optparse::make_option(c("--rds"),
        type = "character", default = NULL,
        help = "Processed rds file of a Seurat object", metavar = "character"
    ),
    optparse::make_option(c("--keywords"),
        type = "character", default = NULL,
        help = "Keywords to slice object", metavar = "character"
    ),
    optparse::make_option(c("--column"),
        type = "character", default = NULL,
        help = "Meta data column name for selecting the slice", metavar = "character"
    ),
    optparse::make_option(c("--output.rds"),
        type = "character", default = "output.rds",
        help = "Output RDS file name [default= %default]", metavar = "character"
    ),
    optparse::make_option(c("--metadata"),
        type = "character", default = NULL,
        help = "Metadata filename", metavar = "character"
    ),
    optparse::make_option(c("--exact"),
        action = "store_true", default = FALSE,
        help = "Exact match, otherwise pattern match"
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



scrna <- readRDS(file = opt$rds)
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

if (!is.null(opt$metadata)) {
    function_read_metadata(opt) -> df
    if (is.null(df)) {
        message("No metadata file provided, skipping metadata merge")
    } else {
        scrna@meta.data %>% dplyr::left_join(df, by = c("orig.ident" = names(df)[1])) -> scrna@meta.data
    }
}




function_subset_by_idents <- function(scrna, opt) {
    if (!is.null(opt$keywords)) {
        keywords <- strsplit(opt$keywords, ",")[[1]]
        patterns <- paste0("(?i)", "(", paste(keywords, collapse = "|"), ")")

        scrna@meta.data %>%
            {
                if (!opt$exact) {
                    if (is.null(opt$column)) dplyr::filter(., if_any(everything(), str_detect, pattern = patterns)) else dplyr::filter(., if_any(opt$column, str_detect, pattern = patterns))
                } else {
                    dplyr::filter(., get(opt$column) %in% keywords)
                }
            } %>%
            rownames() -> cells
        scrna <- subset(scrna, cells = cells)

        return(scrna)
    }
}



function_subset_by_idents(scrna, opt) -> scrna

head(scrna)

saveRDS(scrna, file = opt$output.rds)