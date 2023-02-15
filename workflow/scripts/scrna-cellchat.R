#!/usr/bin/env Rscript


option_list = list(
optparse::make_option(c("--rds"), type="character", default=NULL, 
              help="Processed rds file of a Seurat object", metavar="character"),

optparse::make_option(c("--idents"), type="character", default="seurat_clusters", 
              help="Meta data column name for identities", metavar="character"),

optparse::make_option(c("--species"), type="character", default="human", 
              help="Species", metavar="character"),

optparse::make_option(c("--output"), type="character", default="cellchat.rds", 
              help="Cellchat output file name", metavar="character")


)
 


opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$rds)){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file)", call.=FALSE)
}

require(Seurat)
require(tidyverse)
require(randomcoloR)
require(CellChat)
require(ggalluvial)

if (opt$species == "mouse")
{
	CellChatDB <- CellChatDB.mouse
	ppi=PPI.mouse
}
if (opt$species  == "human")
{
	CellChatDB <- CellChatDB.human
	ppi=PPI.human
}

try({source("workflow/scripts/scrna-functions.R")},silent=TRUE)
try({source(paste0(system("python -c 'import os; import cellsnake; print(os.path.dirname(cellsnake.__file__))'", intern = TRUE),"/scrna/workflow/scripts/scrna-functions.R"))},silent=TRUE)

scrna=readRDS(file = opt$rds)
DefaultAssay(scrna) <- "RNA"



scrna@meta.data %>% dplyr::mutate(identity=paste0("_",as.character(get(opt$idents)))) -> scrna@meta.data


technicalDetailsPathStep2ResolutionCellChat = ""
pathToPlotsStep2ResolutionCellChat = paste0(opt$output.dir,"/")
dir.create(pathToPlotsStep2ResolutionCellChat, showWarnings = F)


cellchat <- createCellChat(object = scrna, group.by = "identity",datatype = "RNA")


cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)

#Part II: Inference of cell-cell communication network

#check if you have data
#cellchat@data
#cellchat@data.signaling
#unique(cellchat@idents)




cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, ppi)
groupSize <- as.numeric(table(cellchat@idents))

cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)

cellchat <- computeCommunProb(cellchat)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)


saveRDS(cellchat,opt$output)