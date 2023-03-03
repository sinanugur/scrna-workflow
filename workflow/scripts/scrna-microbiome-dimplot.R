#!/usr/bin/env Rscript


option_list = list(
        optparse::make_option(c("--rds"), type="character", default=NULL, 
              help="Processed rds file of a Seurat object", metavar="character"),
        optparse::make_option(c("--microbiome.rds"), type="character", default=NULL, 
              help="Microbiome rds file", metavar="character"),
        optparse::make_option(c("--reduction.type"), type="character", default="umap", 
              help="Reduction type, umap or tsne", metavar="character"),
      optparse::make_option(c("--output.plot"), type="character", default="reduction.pdf", 
              help="Plot file name", metavar="character"),
optparse::make_option(c("--taxa"),
    type = "character", default = "genus",
    help = "Taxonomic level", metavar = "character"
  )


)
 


opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$rds)){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file)", call.=FALSE)
}

require(Seurat)
require(tidyverse)
#require(randomcoloR)
require(tidyseurat)
require(viridis)

try({source("workflow/scripts/scrna-functions.R")},silent=TRUE)
try({source(paste0(system("python -c 'import os; import cellsnake; print(os.path.dirname(cellsnake.__file__))'", intern = TRUE),"/scrna/workflow/scripts/scrna-functions.R"))},silent=TRUE)

scrna=readRDS(file = opt$rds)
DefaultAssay(scrna) <- "RNA"


microbiome=readRDS(file = opt$microbiome.rds)



AddMetaData(scrna,microbiome %>% rownames_to_column("barcodes") %>% gather(taxa,umi,-barcodes)  %>% group_by(taxa) %>% mutate(sum=sum(umi,na.rm = T)) %>% ungroup() %>% 
mutate(taxa=ifelse(sum >= min(sort(unique(sum),decreasing = T)[1:11],na.rm=T),paste0(opt$taxa,"_",taxa),paste0(opt$taxa,"_","others"))) %>% 
select(-sum) %>% group_by(barcodes,taxa) %>% summarise(sum=sum(umi,na.rm = T)) %>% ungroup() %>% spread(taxa,sum) %>% column_to_rownames("barcodes")) -> scrna


#p1 <- DimPlot(scrna, reduction = opt$reduction.type, label = TRUE) & theme_cellsnake_classic() & scale_color_manual(values = palette) 

 scrna  %>% select(barcodes=.cell,orig.ident,contains(opt$taxa),starts_with(opt$reduction.type)) %>% gather(taxa,umi,starts_with(opt$taxa)) %>% select(barcodes,orig.ident,x=3,y=4,taxa,umi) %>% replace(is.na(.), 0) %>% ggplot(aes(x=x,y=y,color=log(umi+1))) + geom_point(size=0.2)  +
  labs(color="Log-UMI") + theme(axis.text = element_text(size=12)) + scale_color_viridis(option = "magma",direction = -1,alpha=0.8,na.value="white") + ggthemes::theme_few() + facet_wrap(~taxa) -> p1



ggsave(plot =p1,filename=opt$output.plot,width=9,height=9)

