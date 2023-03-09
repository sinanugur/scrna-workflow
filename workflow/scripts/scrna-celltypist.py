#!/usr/bin/env python

import sys
import scanpy as sc
import pandas as pd
import celltypist



data = sc.read_h5ad(sys.argv[1])
dotplot_file_name=sys.argv[2]
output_analyses_directory=sys.argv[3]
output_xlsx=sys.argv[4]
model=sys.argv[5]
idents= 'seurat_clusters' if sys.argv[6] is None else sys.argv[6]



predictions = celltypist.annotate(data, model = model, majority_voting = True)



adata = predictions.to_adata()
#'seurat_clusters'
dotplot=celltypist.dotplot(predictions, use_as_reference = idents, use_as_prediction = 'majority_voting', show=False,return_fig=True)


dotplot.savefig(dotplot_file_name)

predictions.to_table(output_analyses_directory)

table = pd.crosstab(adata.obs.seurat_clusters, adata.obs.majority_voting)
table.to_excel(output_xlsx)