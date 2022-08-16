#!/usr/bin/env Rscript


packages <- c("tidyverse","optparse","librarian","Seurat","SeuratDisk","patchwork","harmony","DoubletFinder","viridis","clustree","openxlsx","topGO")


installed_packages <- packages %in% rownames(installed.packages())

if (any(installed_packages == FALSE)) {
  
if (!requireNamespace("librarian", quietly = TRUE)) {
  install.packages("librarian") }



librarian::shelf("optparse")
librarian::shelf("tidyverse","patchwork")
librarian::shelf("Seurat")
librarian::shelf('clustree')
librarian::shelf('openxlsx')
librarian::shelf('chris-mcginnis-ucsf/DoubletFinder')
librarian::shelf('viridis')
librarian::shelf('topGO')


if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
librarian::shelf("mojaveazure/seurat-disk") }

}

