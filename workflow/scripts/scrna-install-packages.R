#!/usr/bin/env Rscript

r <- getOption("repos")
r["CRAN"] <- "https://cloud.r-project.org/"
options(repos = r)


packages <- c(
  "tidyverse", "optparse", "librarian", "Seurat", "SeuratDisk", "patchwork",
  "DoubletFinder", "viridis", "clustree", "openxlsx", "topGO", "org.Hs.eg.db", "miQC", "scater", "MultiKParallel", "limma", "ggthemes", "ComplexHeatmap", "CellChat", "NMF", "clusterProfiler",
  "tidyseurat", "SeuratWrappers", "monocle3", "randomcoloR"
)




installed_packages <- packages %in% rownames(installed.packages())



if (any(installed_packages == FALSE)) {
  print("Packages to be installed: ")

  print(packages[!installed_packages])

  if (!requireNamespace("librarian", quietly = TRUE)) {
    install.packages("librarian")
  }

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }

  if (!requireNamespace("Biobase", quietly = TRUE)) {
    BiocManager::install("Biobase")
  }



  librarian::shelf("optparse") # conda installer
  librarian::shelf("tidyverse", "patchwork") # conda installer

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    try({
      remotes::install_version("Seurat", version = "4.3.0")
    })
    librarian::shelf("Seurat")
  }

  librarian::shelf("clustree")
  librarian::shelf("openxlsx")
  librarian::shelf("chris-mcginnis-ucsf/DoubletFinder")
  librarian::shelf("viridis")
  if (!requireNamespace("topGO", quietly = TRUE)) {
    BiocManager::install("topGO")
  }
  librarian::shelf("org.Hs.eg.db")

  if (!requireNamespace("randomcoloR", quietly = TRUE)) {
    install.packages("randomcoloR")
  }

  librarian::shelf("miQC") # conda installer
  # librarian::shelf('scater') #use conda installer
  librarian::shelf("sinanugur/MultiKParallel")
  librarian::shelf("stemangiola/tidyseurat")

  librarian::shelf("limma")
  librarian::shelf("ggthemes")
  librarian::shelf("NMF")
  librarian::shelf("ComplexHeatmap")
  librarian::shelf("clusterProfiler")

  # librarian::shelf('harmony')

  # if (!requireNamespace("cerebroApp", quietly = TRUE)) {
  #  remotes::install_github("romanhaa/cerebroApp")
  # }

  try({
    if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
      remotes::install_github("mojaveazure/seurat-disk")
    }
    if (!requireNamespace("CellChat", quietly = TRUE)) {
      remotes::install_github("sqjin/CellChat")
    }
    if (!requireNamespace("SeuratWrappers", quietly = TRUE)) {
      remotes::install_github("satijalab/seurat-wrappers")
    }
    if (!requireNamespace("monocle3", quietly = TRUE)) {
      remotes::install_github("cole-trapnell-lab/monocle3")
    }
  })
} else {
  print("All packages were installed...OK")
}