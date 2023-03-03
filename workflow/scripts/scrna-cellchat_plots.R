#!/usr/bin/env Rscript


option_list <- list(
  optparse::make_option(c("--rds"),
    type = "character", default = NULL,
    help = "Cellchat object", metavar = "character"
  ),
  optparse::make_option(c("--output.dir"),
    type = "character", default = NULL,
    help = "Output directory", metavar = "character"
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
require(CellChat)
require(ggalluvial)

cellchat <- readRDS(file = opt$rds)




technicalDetailsPathStep2ResolutionCellChat <- ""
pathToPlotsStep2ResolutionCellChat <- paste0(opt$output.dir, "/")
dir.create(pathToPlotsStep2ResolutionCellChat, showWarnings = F)



groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1, 2), xpd = TRUE)
fName <- "Number-of-interactions"
pdf(paste(pathToPlotsStep2ResolutionCellChat, fName, ".pdf", sep = ""), 10, 10)
print(netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Number of interactions"))
dev.off()

fName <- "Interaction-strength"
pdf(paste(pathToPlotsStep2ResolutionCellChat, fName, ".pdf", sep = ""), 10, 10)
print(netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Interaction weights/strength"))
dev.off()

pathToPlotsStep2ResolutionCellChatClusters <- paste0(pathToPlotsStep2ResolutionCellChat, "1-clusters", "/")
dir.create(pathToPlotsStep2ResolutionCellChatClusters, showWarnings = F)


mat <- cellchat@net$weight
par(mfrow = c(3, 4), xpd = TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  fName <- stringr::str_replace_all(paste0("edge-weights-for-", rownames(mat)[i]), "/", "_")
  pdf(paste(pathToPlotsStep2ResolutionCellChatClusters, fName, ".pdf", sep = ""), 10, 10)
  print(netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i]))
  dev.off()
}

pathToPlotsStep2ResolutionCellChatSignPath <- paste0(pathToPlotsStep2ResolutionCellChat, "2-significant-pathways", "/")
dir.create(pathToPlotsStep2ResolutionCellChatSignPath, showWarnings = F)

# All the signaling pathways showing significant communications
allPathways <- cellchat@netP$pathways
for (pathways.show in allPathways)
{
  pathToPlotsStep2ResolutionCellChatSignPathCicrcle <- paste0(pathToPlotsStep2ResolutionCellChatSignPath, "circle-plots", "/")
  dir.create(pathToPlotsStep2ResolutionCellChatSignPathCicrcle, showWarnings = F)

  # Hierarchy plot
  # Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells
  vertex.receiver <- seq(1, 4) # a numeric vector.
  netVisual_aggregate(cellchat, signaling = pathways.show, vertex.receiver = vertex.receiver)
  # Circle plot
  par(mfrow = c(1, 1))
  fName <- pathways.show
  pdf(paste(pathToPlotsStep2ResolutionCellChatSignPathCicrcle, fName, ".pdf", sep = ""), 10, 10)
  print(netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle"))
  dev.off()

  pathToPlotsStep2ResolutionCellChatSignPathchord <- paste0(pathToPlotsStep2ResolutionCellChatSignPath, "chord-plots", "/")
  dir.create(pathToPlotsStep2ResolutionCellChatSignPathchord, showWarnings = F)
  # Chord diagram
  par(mfrow = c(1, 1))
  pdf(paste(pathToPlotsStep2ResolutionCellChatSignPathchord, fName, ".pdf", sep = ""), 10, 10)
  print(netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord"))
  dev.off()

  pathToPlotsStep2ResolutionCellChatSignPathHeatmap <- paste0(pathToPlotsStep2ResolutionCellChatSignPath, "heatmap-plots", "/")
  dir.create(pathToPlotsStep2ResolutionCellChatSignPathHeatmap, showWarnings = F)
  # Heatmap
  par(mfrow = c(1, 1))

  pdf(paste(pathToPlotsStep2ResolutionCellChatSignPathHeatmap, fName, ".pdf", sep = ""), 10, 10)
  print(netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds"))
  dev.off()
}


#> Do heatmap based on a single object


# Chord diagram
# group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
# names(group.cellType) <- levels(cellchat@idents)
# netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
#> Plot the aggregated cell-cell communication network at the signaling pathway level


pathToPlotsStep2ResolutionCellChatContributionOfLRpairs <- paste0(pathToPlotsStep2ResolutionCellChat, "3-ligand-receptor-pairs", "/")
dir.create(pathToPlotsStep2ResolutionCellChatContributionOfLRpairs, showWarnings = F)
fName <- "all-pathways-contribution"
tryCatch(
  {
    pdf(paste(pathToPlotsStep2ResolutionCellChatContributionOfLRpairs, fName, ".pdf", sep = ""), 5, 20)
    print(netAnalysis_contribution(cellchat, signaling = allPathways))
    dev.off()
  },
  error = function(e) {
    dev.off()
  }
)

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
for (pathways.show in allPathways)
{
  pathToPlotsStep2ResolutionCellChatContributionOfLRpairsCircle <- paste0(pathToPlotsStep2ResolutionCellChatContributionOfLRpairs, "circle", "/")
  dir.create(pathToPlotsStep2ResolutionCellChatContributionOfLRpairsCircle, showWarnings = F)

  pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
  LR.show <- pairLR.CXCL[1, ] # show one ligand-receptor pair
  # Hierarchy plot
  vertex.receiver <- seq(1, 4) # a numeric vector
  fName <- pathways.show
  pdf(paste(pathToPlotsStep2ResolutionCellChatContributionOfLRpairsCircle, fName, ".pdf", sep = ""), 10, 10)
  print(netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle"))
  dev.off()

  pathToPlotsStep2ResolutionCellChatContributionOfLRpairsChord <- paste0(pathToPlotsStep2ResolutionCellChatContributionOfLRpairs, "chord", "/")
  dir.create(pathToPlotsStep2ResolutionCellChatContributionOfLRpairsChord, showWarnings = F)

  #> [[1]]
  # Chord diagram
  pdf(paste(pathToPlotsStep2ResolutionCellChatContributionOfLRpairsChord, fName, ".pdf", sep = ""), 10, 10)
  print(netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord"))
  dev.off()
}

pathToPlotsStep2ResolutionCellChatCellCellC <- paste0(pathToPlotsStep2ResolutionCellChat, "4-cell-cell-communication", "/")
dir.create(pathToPlotsStep2ResolutionCellChatCellCellC, showWarnings = F)

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
fName <- "significant-interactions-for-cell-groups"
pdf(paste(pathToPlotsStep2ResolutionCellChatCellCellC, fName, ".pdf", sep = ""), 50, 20)
print(netVisual_bubble(cellchat, remove.isolate = FALSE))
dev.off()
# netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)


# show all the significant interactions (L-R pairs) associated with certain signaling pathways
fName <- "significant-interactoions-for-singaling-pathways"
pdf(paste(pathToPlotsStep2ResolutionCellChatCellCellC, fName, ".pdf", sep = ""), 50, 20)
print(netVisual_bubble(cellchat, signaling = allPathways, remove.isolate = FALSE))
dev.off()

# show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
# pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","FGF"))
# netVisual_bubble(cellchat, sources.use = c(3,4), targets.use = c(5:8), pairLR.use = pairLR.use, remove.isolate = TRUE)
#> Comparing communications on a single object
#>

pathToPlotsStep2ResolutionCellChatCellCellCAllInteractions <- paste0(pathToPlotsStep2ResolutionCellChatCellCellC, "interactions", "/")
dir.create(pathToPlotsStep2ResolutionCellChatCellCellCAllInteractions, showWarnings = F)

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB

if (FALSE) { # commented out because it takes a lot of time
  for (g1 in unique(cellchat@meta$identity))
  {
    for (g2 in unique(cellchat@meta$identity))
    {
      tryCatch(
        {
          fName <- paste0("from-", g1, "-to-", g2)
          pdf(paste(pathToPlotsStep2ResolutionCellChatCellCellCAllInteractions, fName, ".pdf", sep = ""), 50, 20)
          print(netVisual_chord_gene(cellchat, sources.use = g1, targets.use = g2, lab.cex = 0.5, legend.pos.y = 30))
          dev.off()
        },
        error = function(cond) {
          dev.off()
        },
        finally = {}
      )
    }
  }

  pathToPlotsStep2ResolutionCellChatCellCellCAllInteractionsWithP <- paste0(pathToPlotsStep2ResolutionCellChatCellCellC, "interactions-with-pathways", "/")
  dir.create(pathToPlotsStep2ResolutionCellChatCellCellCAllInteractionsWithP, showWarnings = F)

  for (g1 in unique(cellchat@meta$identity))
  {
    for (g2 in unique(cellchat@meta$identity))
    {
      for (pathways.show in allPathways)
      {
        tryCatch(
          {
            fName <- paste0("from-", g1, "-to-", g2, "-pathway-", pathways.show)
            pdf(paste(pathToPlotsStep2ResolutionCellChatCellCellCAllInteractionsWithP, fName, ".pdf", sep = ""), 50, 20)
            print(netVisual_chord_gene(cellchat, sources.use = g1, targets.use = g2, signaling = pathways.show, lab.cex = 0.5, legend.pos.y = 30))
            dev.off()
          },
          error = function(cond) {
            dev.off()
          },
          finally = {}
        )
      }
    }
  }
}

pathToPlotsStep2ResolutionGE <- paste0(pathToPlotsStep2ResolutionCellChat, "5-gene-expression", "/")
dir.create(pathToPlotsStep2ResolutionGE, showWarnings = F)

for (pathways.show in allPathways)
{
  fName <- paste0("gene-expression", "-pathway-", pathways.show)
  pdf(paste(pathToPlotsStep2ResolutionGE, fName, ".pdf", sep = ""), 20, 10)
  print(plotGeneExpression(cellchat, signaling = pathways.show))
  dev.off()
}

pathToPlotsStep2ResolutionCC <- paste0(pathToPlotsStep2ResolutionCellChat, "6-system-analysis-of-cc-communication-network", "/")
dir.create(pathToPlotsStep2ResolutionCC, showWarnings = F)

pathToPlotsStep2ResolutionCCNCS <- paste0(pathToPlotsStep2ResolutionCC, "network-centrality-scores/")
dir.create(pathToPlotsStep2ResolutionCCNCS, showWarnings = F)
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
for (pathways.show in allPathways)
{
  fName <- pathways.show
  pdf(paste(pathToPlotsStep2ResolutionCCNCS, fName, ".pdf", sep = ""), 5, 5)
  print(netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, font.size = 10))
  dev.off()
}

pathToPlotsStep2ResolutionCCDSR <- paste0(pathToPlotsStep2ResolutionCC, "dominant-senders-and-receivers/")
dir.create(pathToPlotsStep2ResolutionCCDSR, showWarnings = F)

# Visualize the dominant senders (sources) and receivers (targets) in a 2D space

# We also provide another intutive way to visualize the dominant senders (sources) and receivers (targets) in a 2D space using scatter plot.

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
fName <- "aggregated-results"
pdf(paste(pathToPlotsStep2ResolutionCC, fName, ".pdf", sep = ""), 10, 10)
print(netAnalysis_signalingRole_scatter(cellchat))
dev.off()

pathToPlotsStep2ResolutionCCDSRSP <- paste0(pathToPlotsStep2ResolutionCCDSR, "scatter/")
dir.create(pathToPlotsStep2ResolutionCCDSRSP, showWarnings = F)

pathToPlotsStep2ResolutionCCDH <- paste0(pathToPlotsStep2ResolutionCCDSR, "heatmap/")
dir.create(pathToPlotsStep2ResolutionCCDH, showWarnings = F)

for (pathways.show in allPathways)
{
  #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  # Signaling role analysis on the cell-cell communication networks of interest

  fName <- pathways.show

  tryCatch(
    {
      pdf(paste(pathToPlotsStep2ResolutionCCDSRSP, fName, ".pdf", sep = ""), 10, 20)
      print(netAnalysis_signalingRole_scatter(cellchat, signaling = pathways.show))
      dev.off()
    },
    error = function(cond) {
      dev.off()
    },
    finally = {}
  )

  tryCatch(
    {
      # Signaling role analysis on the cell-cell communication networks of interest
      pdf(paste(pathToPlotsStep2ResolutionCCDH, fName, ".pdf", sep = ""), 10, 20)
      print(netAnalysis_signalingRole_heatmap(cellchat, signaling = pathways.show))
      dev.off()
    },
    error = function(cond) {
      dev.off()
    },
    finally = {}
  )
}

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
fName <- "aggregated-results-incoming-outgoing"
pdf(paste(pathToPlotsStep2ResolutionCC, fName, ".pdf", sep = ""), 10, 5)
print(ht1 + ht2)
dev.off()



# Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together

pathToPlotsStep2ResolutionOCPWSC <- paste0(pathToPlotsStep2ResolutionCellChat, "7-communication-pattern-of-secreting-cells/")
dir.create(pathToPlotsStep2ResolutionOCPWSC, showWarnings = F)

# write("In outgoing: You look at Cophenetic and Silhouette values and when they begin to drop suddenly then you look at the number (one before dropping), opposite with icoming pattern", paste0(pathToPlotsStep2ResolutionOCPWSC, "instruction.txt"))

library(NMF)
library(ggalluvial)

pathToPlotsStep2ResolutionOCPWSCOut <- paste0(pathToPlotsStep2ResolutionOCPWSC, "outgoing/")
dir.create(pathToPlotsStep2ResolutionOCPWSCOut, showWarnings = F)

out <- selectK(cellchat, pattern = "outgoing")

fName <- "outgoing-communication-pattern"
pdf(paste(pathToPlotsStep2ResolutionOCPWSCOut, fName, ".pdf", sep = ""), 10, 5)
print(out)
dev.off()

pathToPlotsStep2ResolutionOCPWSCHeatmap <- paste0(pathToPlotsStep2ResolutionOCPWSCOut, "heatmap/")
dir.create(pathToPlotsStep2ResolutionOCPWSCHeatmap, showWarnings = F)

pathToPlotsStep2ResolutionOCPWSCRiver <- paste0(pathToPlotsStep2ResolutionOCPWSCOut, "river/")
dir.create(pathToPlotsStep2ResolutionOCPWSCRiver, showWarnings = F)

pathToPlotsStep2ResolutionOCPWSCDot <- paste0(pathToPlotsStep2ResolutionOCPWSCOut, "dot/")
dir.create(pathToPlotsStep2ResolutionOCPWSCDot, showWarnings = F)


for (nPatterns in c(1:10))
{
  tryCatch(
    {
      fName <- paste0("pattern-", nPatterns)
      cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
      pdf(paste(pathToPlotsStep2ResolutionOCPWSCHeatmap, fName, ".pdf", sep = ""), 10, 10)
      print(identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns))
      dev.off()

      pdf(paste(pathToPlotsStep2ResolutionOCPWSCRiver, fName, ".pdf", sep = ""), 10, 10)
      print(netAnalysis_river(cellchat, pattern = "outgoing"))
      dev.off()

      pdf(paste(pathToPlotsStep2ResolutionOCPWSCDot, fName, ".pdf", sep = ""), 10, 10)
      print(netAnalysis_dot(cellchat, pattern = "outgoing"))
      dev.off()
    },
    error = function(cond) {},
    finally = {}
  )
}




pathToPlotsStep2ResolutionOCPWSCIn <- paste0(pathToPlotsStep2ResolutionOCPWSC, "incoming/")
dir.create(pathToPlotsStep2ResolutionOCPWSCIn, showWarnings = F)

out <- selectK(cellchat, pattern = "incoming")

fName <- "incoming-communication-pattern"
pdf(paste(pathToPlotsStep2ResolutionOCPWSCIn, fName, ".pdf", sep = ""), 10, 10)
print(out)
dev.off()

pathToPlotsStep2ResolutionOCPWSCHeatmap <- paste0(pathToPlotsStep2ResolutionOCPWSCIn, "heatmap/")
dir.create(pathToPlotsStep2ResolutionOCPWSCHeatmap, showWarnings = F)


pathToPlotsStep2ResolutionOCPWSCRiver <- paste0(pathToPlotsStep2ResolutionOCPWSCIn, "river/")
dir.create(pathToPlotsStep2ResolutionOCPWSCRiver, showWarnings = F)

pathToPlotsStep2ResolutionOCPWSCDot <- paste0(pathToPlotsStep2ResolutionOCPWSCIn, "dot/")
dir.create(pathToPlotsStep2ResolutionOCPWSCDot, showWarnings = F)


for (nPatterns in 1:10)
{
  tryCatch(
    {
      fName <- paste0("pattern-", nPatterns)
      cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
      pdf(paste(pathToPlotsStep2ResolutionOCPWSCHeatmap, fName, ".pdf", sep = ""), 10, 10)
      print(identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns))
      dev.off()

      pdf(paste(pathToPlotsStep2ResolutionOCPWSCRiver, fName, ".pdf", sep = ""), 10, 10)
      print(netAnalysis_river(cellchat, pattern = "incoming"))
      dev.off()

      pdf(paste(pathToPlotsStep2ResolutionOCPWSCDot, fName, ".pdf", sep = ""), 10, 10)
      print(netAnalysis_dot(cellchat, pattern = "incoming"))
      dev.off()
    },
    error = function(cond) {},
    finally = {}
  )
}