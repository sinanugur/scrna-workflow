#!/usr/bin/env Rscript
option_list = list(

    optparse::make_option(c("--xlsx"), type="character", default=NULL, 
              help="Excel table of markers", metavar="character"),
    optparse::make_option(c("--output"), type="character", default=NULL, 
              help="Output excel file name", metavar="character"),
    optparse::make_option(c("--ontology"), type="character", default="BP", 
              help="GO ontology, possible values BP, CC or MF", metavar="character"),
    optparse::make_option(c("--algorithm"), type="character", default="weight01", 
              help="Algorithm", metavar="character"),
        optparse::make_option(c("--mapping"), type="character", default="org.Hs.eg.db", 
              help="Mapping", metavar="character"),
                 optparse::make_option(c("--statistics"), type="character", default="ks", 
              help="Statistics", metavar="character")


)
 


opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$xlsx) || is.null(opt$output) ){
  optparse::print_help(opt_parser)
  stop("Arguments must be supplied", call.=FALSE)
}



require(Seurat)
require(tidyverse)
require(topGO)
#source("workflow/scripts/scrna-functions.R")

function_gofunc=function(df,algorithm="weight01",statistics="ks",mapping="org.Hs.eg.db",ID="symbol",ontology = "BP") {
  

geneList <- df$p_val
names(geneList) <- df$gene
# Create topGOData object
GOdata <- new("topGOdata",
              ontology = ontology,
              allGenes = geneList,
              geneSelectionFun = function(x)x,
              annot = annFUN.org, mapping = mapping,ID=ID,nodeSize = 10)

resultsKS <- runTest(GOdata,algorithm = algorithm,statistic = statistics)


tab <- GenTable(GOdata, raw.p.value = resultsKS, topNodes = length(resultsKS@score), numChar = 120)

return(tab)
}



All_Features=openxlsx::read.xlsx(opt$xlsx)


All_Features %>% group_by(cluster) %>% split(.$cluster) %>% map(~function_gofunc(.)) %>% bind_rows(.id="cluster") -> go_enrichment

openxlsx::write.xlsx(go_enrichment,file=opt$output)






