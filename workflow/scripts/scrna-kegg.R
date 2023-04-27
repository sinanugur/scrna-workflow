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
  optparse::make_option(c("--output.kegg"),
    type = "character", default = NULL,
    help = "Output kegg excel file name", metavar = "character"
  ),
  optparse::make_option(c("--output.mkegg"),
    type = "character", default = NULL,
    help = "Output kegg excel file name", metavar = "character"
  ),
  optparse::make_option(c("--output.gse"),
    type = "character", default = NULL,
    help = "Output gse excel file name", metavar = "character"
  ),
  optparse::make_option(c("--output.mgse"),
    type = "character", default = NULL,
    help = "Output gse excel file name", metavar = "character"
  ),
  optparse::make_option(c("--mapping"),
    type = "character", default = "org.Hs.eg.db",
    help = "Mapping", metavar = "character"
  ),
  optparse::make_option(c("--organism"),
    type = "character", default = "hsa",
    help = "Organism code, look at https://www.genome.jp/kegg/catalog/org_list.html", metavar = "character"
  ),
  optparse::make_option(c("--pval"),
    type = "double", default = 0.05,
    help = "P value treshold [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--logfc.treshold"),
    type = "double", default = 1,
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


function_enrichment_kegg_singlecell <- function(results, p = 0.05, f = 1.5) {
  print(results %>% distinct(cluster) %>% dplyr::pull())
  results %>%
    as.data.frame() %>%
    dplyr::filter(p_val_adj < p) %>%
    arrange(desc(avg_log2FC)) %>%
    dplyr::select(gene, avg_log2FC) %>%
    dplyr::mutate(GeneID = mapIds(get(opt$mapping), keys = gene, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")) %>%
    dplyr::filter(!is.na(GeneID), !is.na(avg_log2FC), !duplicated(GeneID)) %>%
    dplyr::select(3, 2) %>%
    deframe() -> geneList


  gene <- names(geneList)[abs(geneList) > f] # as recommended by the documentation

  tryCatch(
    {
      kk1 <- clusterProfiler::enrichKEGG(
        gene = gene,
        organism = opt$organism,
        pAdjustMethod = "fdr",
        minGSSize = 2,
        pvalueCutoff = 1
      )
    },
    error = function(e) {
      kk1 <- NULL
    }
  ) -> kk1

  tryCatch(
    {
      kk2 <- clusterProfiler::gseKEGG(
        geneList = geneList,
        organism = opt$organism,
        pvalueCutoff = 1,
        pAdjustMethod = "fdr",
        minGSSize = 2,
        eps = 0,
        verbose = FALSE
      )
    },
    error = function(e) {
      kk2 <- NULL
    }
  ) -> kk2

  tryCatch(
    {
      kk3 <- clusterProfiler::enrichMKEGG(
        gene = gene,
        organism = opt$organism,
        pvalueCutoff = 1,
        minGSSize = 2,
        pAdjustMethod = "fdr",
      )
    },
    error = function(e) {
      kk3 <- NULL
    }
  ) -> kk3

  tryCatch(
    {
      kk4 <- clusterProfiler::gseMKEGG(
        geneList = geneList,
        organism = opt$organism,
        minGSSize = 2,
        keyType = "kegg",
        pAdjustMethod = "fdr",
        pvalueCutoff = 1
      )
    },
    error = function(e) {
      kk4 <- NULL
    }
  ) -> kk4

  return(list(kk1, kk2, kk3, kk4))
}


All_Features %>%
  split(.$cluster) %>%
  purrr::map(~ function_enrichment_kegg_singlecell(., p = opt$pval, f = opt$logfc.treshold)) -> all_kegg_results



saveRDS(all_kegg_results, opt$output.rds)

openxlsx::write.xlsx(all_kegg_results %>% keep(~ !is.null(.[[1]])) %>% map(~ .[[1]]@result) %>% bind_rows(.id = "cluster") %>% as_tibble(), opt$output.kegg)
openxlsx::write.xlsx(all_kegg_results %>% keep(~ !is.null(.[[2]])) %>% map(~ .[[2]]@result) %>% bind_rows(.id = "cluster") %>% as_tibble(), opt$output.gse)
openxlsx::write.xlsx(all_kegg_results %>% keep(~ !is.null(.[[3]])) %>% map(~ .[[3]]@result) %>% bind_rows(.id = "cluster") %>% as_tibble(), opt$output.mkegg)
openxlsx::write.xlsx(all_kegg_results %>% keep(~ !is.null(.[[4]])) %>% map(~ .[[4]]@result) %>% bind_rows(.id = "cluster") %>% as_tibble(), opt$output.mgse)