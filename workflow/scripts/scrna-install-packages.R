#!/usr/bin/env Rscript


packages <- c("tidyverse","optparse","librarian","Seurat","SeuratDisk","patchwork","harmony")


installed_packages <- packages %in% rownames(installed.packages())

if (any(installed_packages == FALSE)) {
  
if (!requireNamespace("librarian", quietly = TRUE)) {
  install.packages("librarian") }



librarian::shelf("optparse")
librarian::shelf("tidyverse","patchwork")
librarian::shelf("Seurat")



if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
  librarian::shelf("mojaveazure/seurat-disk") }


}

