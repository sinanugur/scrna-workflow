#!/usr/bin/env Rscript
option_list <- list(
  optparse::make_option(c("--xlsx"),
    type = "character", default = NULL,
    help = "Excel table of markers", metavar = "character"
  ),
  optparse::make_option(c("--output.rds"),
    type = "character", default = NULL,
    help = "Output RDS file name", metavar = "character"
  ),
  optparse::make_option(c("--output.go"),
    type = "character", default = NULL,
    help = "Output kegg excel file name", metavar = "character"
  ),
  optparse::make_option(c("--output.gse"),
    type = "character", default = NULL,
    help = "Output gse excel file name", metavar = "character"
  ),
  optparse::make_option(c("--mapping"),
    type = "character", default = "org.Hs.eg.db",
    help = "Mapping", metavar = "character"
  ),
  optparse::make_option(c("--pval"),
    type = "double", default = 0.05,
    help = "P value treshold [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--logfc.treshold"),
    type = "double", default = 1.5,
    help = "LogFC [default= %default]", metavar = "character"
  )
)
# require(clusterProfiler)
require(tidyverse)

try({
  if (!requireNamespace(opt$mapping, quietly = TRUE)) {
    BiocManager::install(opt$mapping)
  }
})


opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

require(opt$mapping, character.only = T)

if (is.null(opt$xlsx)) {
  optparse::print_help(opt_parser)
  stop("Arguments must be supplied", call. = FALSE)
}


All_Features <- openxlsx::read.xlsx(opt$xlsx)


function_enrichment_go_singlecell <- function(results, p = 0.05, f = 1.5) {
  print(results %>% distinct(cluster) %>% pull())
  results %>%
    as.data.frame() %>%
    dplyr::filter(p_val_adj < p) %>%
    arrange(desc(avg_log2FC)) %>%
    dplyr::select(gene, avg_log2FC) %>%
    dplyr::mutate(GeneID = mapIds(get(opt$mapping), keys = gene, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")) %>%
    dplyr::filter(!is.na(GeneID), !is.na(avg_log2FC), !duplicated(GeneID)) %>%
    dplyr::select(3, 2) %>%
    deframe() -> geneList


  gene <- names(geneList)[abs(geneList) > f]

  tryCatch(
    {
      kk <- clusterProfiler::enrichGO(
        gene = gene,
        universe = names(geneList),
        OrgDb = get(opt$mapping),
        ont = "ALL",
        pAdjustMethod = "fdr",
        pvalueCutoff = 0.05,
        readable = TRUE
      )
    },
    error = function(e) {
      kk <- NULL
    }
  ) -> kk

  tryCatch(
    {
      kk2 <- clusterProfiler::gseGO(
        geneList = names(geneList),
        OrgDb = get(opt$mapping),
        ont = "ALL",
        minGSSize = 2,
        maxGSSize = 500,
        pvalueCutoff = 0.05,
        verbose = FALSE
      )
    },
    error = function(e) {
      kk2 <- NULL
    }
  ) -> kk2





  return(list(kk, kk2))
}


All_Features %>%
  split(.$cluster) %>%
  purrr::map(~ function_enrichment_go_singlecell(., p = opt$pval, f = opt$logfc.treshold)) -> all_go_results



saveRDS(all_go_results, opt$output.rds)

openxlsx::write.xlsx(all_go_results %>% keep(~ !is.null(.[[1]])) %>% map(~ .[[1]]@result) %>% bind_rows(.id = "cluster") %>% as_tibble(), opt$output.go)
openxlsx::write.xlsx(all_go_results %>% keep(~ !is.null(.[[2]])) %>% map(~ .[[2]]@result) %>% bind_rows(.id = "cluster") %>% as_tibble(), opt$output.gse)