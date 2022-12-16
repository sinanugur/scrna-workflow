#!/usr/bin/env Rscript
option_list <- list(
  optparse::make_option(c("--xlsx"),
    type = "character", default = NULL,
    help = "Excel table of markers", metavar = "character"
  ),
  optparse::make_option(c("--output"),
    type = "character", default = NULL,
    help = "Output excel file name", metavar = "character"
  ),
  optparse::make_option(c("--ontology"),
    type = "character", default = "BP",
    help = "GO ontology, possible values BP, CC or MF", metavar = "character"
  ),
  optparse::make_option(c("--algorithm"),
    type = "character", default = "weight01",
    help = "Algorithm", metavar = "character"
  ),
  optparse::make_option(c("--mapping"),
    type = "character", default = "org.Hs.eg.db",
    help = "Mapping", metavar = "character"
  ),
  optparse::make_option(c("--statistics"),
    type = "character", default = "ks",
    help = "Statistics", metavar = "character"
  )
)



opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

if (is.null(opt$xlsx) || is.null(opt$output)) {
  optparse::print_help(opt_parser)
  stop("Arguments must be supplied", call. = FALSE)
}


All_Features <- openxlsx::read.xlsx(opt$xlsx)

res <- All_Features
original_gene_list <- res$avg_log2FC
names(original_gene_list) <- res$gene
gene_list <- na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list <- sort(gene_list, decreasing = TRUE)

# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids <- bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# remove duplicate IDS (here we use "SYMBOL", but it should be whatever was selected as keyType)
dedup_ids <- ids[!duplicated(ids[c("SYMBOL")]), ]

# Create a new dataframe df2 which has the respective entrez IDs for the gene symbols
colnames(dedup_ids) <- c("gene", "EntrezID")
df2 <- merge(res, dedup_ids, by = "gene")

# Create a vector of the gene universe
kegg_gene_list <- df2$avg_log2FC
names(kegg_gene_list) <- df2$EntrezID

kegg_gene_list <- na.omit(kegg_gene_list)

# sort the list in decreasing order
kegg_gene_list <- sort(kegg_gene_list, decreasing = TRUE)

kegg_organism <- "hsa"

kk2 <- gseKEGG(
  geneList = kegg_gene_list,
  organism = kegg_organism,
  minGSSize = 3,
  maxGSSize = 800,
  pvalueCutoff = 0.05,
  pAdjustMethod = "none",
  keyType = "ncbi-geneid"
)

saveRDS(kk2, file = "kk2.rds")
# out = as.matrix(kk2@result)
# return(out)