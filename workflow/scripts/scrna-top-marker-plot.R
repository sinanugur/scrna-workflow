#!/usr/bin/env Rscript
option_list <- list(
    optparse::make_option(c("--xlsx"),
        type = "character", default = NULL,
        help = "Excel table of markers for input", metavar = "character"
    ),
    optparse::make_option(c("--output.plot"),
        type = "character", default = "output.pdf",
        help = "Output plot directory", metavar = "character"
    )
)



opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# try({source("workflow/scripts/scrna-functions.R")},silent=TRUE)
# try({source(paste0(system("python -c 'import os; import cellsnake; print(os.path.dirname(cellsnake.__file__))'", intern = TRUE),"/scrna/workflow/scripts/scrna-functions.R"))},silent=TRUE)

if (is.null(opt$xlsx)) {
    optparse::print_help(opt_parser)
    stop("At least one argument must be supplied (rds file)", call. = FALSE)
}

require(tidyverse)



Positive_Features <- openxlsx::read.xlsx(opt$xlsx)


Positive_Features %>%
    group_by(cluster) %>%
    slice_max(n = 20, order_by = avg_log2FC) -> df


Positive_Features %>%
    distinct(cluster) %>%
    pull() -> clusters



pdf(opt$output.plot, width = 6, height = 6)
for (i in clusters) {
    df %>% dplyr::filter(cluster %in% i) -> df2
    maxFC <- (df2 %>% pull(avg_log2FC) %>% max()) + 1

    try({
        df2 %>%
            dplyr::mutate(n = dense_rank(dplyr::desc(avg_log2FC))) %>%
            ggplot(aes(x = n, y = avg_log2FC, label = gene)) +
            geom_text(angle = 75, size = 4) +
            ggthemes::theme_few() +
            ylim(c(0, maxFC)) +
            coord_cartesian(clip = "off", expand = TRUE) +
            ggtitle("Top markers", subtitle = paste(i, "vs all")) -> p1
        print(p1)
    })
}
dev.off()