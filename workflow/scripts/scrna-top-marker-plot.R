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

maxFC <- (df %>% pull(avg_log2FC) %>% max()) + 1

df %>%
    group_by(cluster) %>%
    dplyr::mutate(n = dense_rank(dplyr::desc(avg_log2FC))) %>%
    ggplot(aes(x = n, y = avg_log2FC, label = gene)) +
    geom_text(angle = 75, size = 3) +
    facet_wrap(~cluster, ncol = 4) +
    # ggthemes::theme_few() +
    theme(strip.text = element_text(size = 12)) +
    ylim(c(0, maxFC)) +
    coord_cartesian(clip = "off", expand = TRUE) +
    ggtitle("Top markers") -> p1

Positive_Features %>%
    distinct(cluster) %>%
    pull() %>%
    length() -> n

ggsave(opt$output.plot, p1, height = (2.1 + ceiling(n / 4) * 2.1), width = 9.5, useDingbats = TRUE)