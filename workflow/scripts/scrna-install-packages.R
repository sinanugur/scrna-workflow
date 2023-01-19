#!/usr/bin/env Rscript

r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)


packages <- c("tidyverse","optparse","librarian","Seurat","SeuratDisk","patchwork","harmony","DoubletFinder","viridis","clustree","openxlsx","topGO","org.Hs.eg.db","cerebroApp","miQC","scater")


installed_packages <- packages %in% rownames(installed.packages())

if (any(installed_packages == FALSE)) {
  
if (!requireNamespace("librarian", quietly = TRUE)) {
  install.packages("librarian") }

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager") }

if (!requireNamespace("Biobase", quietly = TRUE)) {
  BiocManager::install("Biobase") }



librarian::shelf("optparse")
librarian::shelf("tidyverse","patchwork")
librarian::shelf("Seurat")
librarian::shelf('clustree')
librarian::shelf('openxlsx')
librarian::shelf('chris-mcginnis-ucsf/DoubletFinder')
librarian::shelf('viridis')
librarian::shelf('topGO')
librarian::shelf('org.Hs.eg.db')
librarian::shelf('randomcoloR')
librarian::shelf('miQC')
librarian::shelf('scater')


if (!requireNamespace("cerebroApp", quietly = TRUE)) {
ibrarian::shelf('romanhaa/cerebroApp') #will check later
}

if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
librarian::shelf("mojaveazure/seurat-disk") }

}

