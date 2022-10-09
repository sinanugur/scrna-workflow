from collections import defaultdict
from yaml import load
import os


def input_function(wildcards):
    if os.path.isfile("data/" + wildcards.sample + "/filtered_feature_bc_matrix/matrix.mtx.gz"):
        return("data/" + wildcards.sample + "/filtered_feature_bc_matrix/")
    elif os.path.isfile("data/" + wildcards.sample + "/filtered_feature_bc_matrix.h5"):
        return("data/" + wildcards.sample + "/")
    elif os.path.isfile("data/" + wildcards.sample + "/raw_feature_bc_matrix/matrix.mtx.gz"):
        return("data/" + wildcards.sample + "/raw_feature_bc_matrix/")
    else:
        return("data/" + wildcards.sample + "/outs/raw_feature_bc_matrix/")


rule rds_params:
    input:
        #"data/{sample}/raw_feature_bc_matrix/"
        input_function
    output:
        rds="analyses/raw/{sample}/" + f"{paramspace.wildcard_pattern}.rds",
        before="results/{sample}/" + f"{paramspace.wildcard_pattern}" + "/technicals/before-qc-trimming-violinplot.pdf",
        after="results/{sample}/" + f"{paramspace.wildcard_pattern}" + "/technicals/after-qc-trimming-violinplot.pdf"
    params:
        paramaters=paramspace.instance,
    shell:
        """workflow/scripts/scrna-read-qc.R --data.dir {input} --output.rds {output.rds} --sampleid {wildcards.sample} --percent.mt {params.paramaters[MT]} --min.features {min_features} --min.cells {min_cells} --before.violin.plot {output.before} --after.violin.plot {output.after}"""

"""
rule rds:
    input:
        #"data/{sample}/raw_feature_bc_matrix/"
        input_function
    output:
        "analyses/raw/{sample}.rds"

    shell:
        "workflow/scripts/scrna-read-qc.R --data.dir {input} --sampleid {wildcards.sample} --percent.mt {percent_mt} --min.features {min_features} --min.cells {min_cells}"
"""

rule clustree:
    input:
        "analyses/raw/{sample}/" + f"{paramspace.wildcard_pattern}.rds"
    output:
        clustree="results/{sample}/" + f"{paramspace.wildcard_pattern}" + "/clusteringTree/clusteringTree.pdf",
        heatmap="results/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/technicals/DimHeatMap_plot.pdf",
        hvfplot="results/{sample}/" + f"{paramspace.wildcard_pattern}" + "/technicals/highly-variable-features.pdf",
        jackandelbow="results/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/technicals/JackandElbow_plot.pdf"
    shell:
        "workflow/scripts/scrna-clusteringtree.R --rds {input} --output {output.clustree} --heatmap {output.heatmap} --hvfplot {output.hvfplot} --jackandelbow {output.jackandelbow}"



rule normalization_pca_rds:
    input:
        "analyses/raw/{sample}/" + f"{paramspace.wildcard_pattern}.rds"
    output:
        rds="analyses/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds",
        xlsx="results/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/number-of-cells-per-cluster.xlsx",
        pca="results/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/pca.plot.pdf"
    params:
        paramaters=paramspace.instance,
        doublet_filter=doublet_filter,
        umap_plot=umap_plot,
        tsne_plot=tsne_plot
    shell:
        "workflow/scripts/scrna-normalization-pca.R --rds {input} {params.doublet_filter} --normalization.method {normalization_method} "
        "--scale.factor {scale_factor} --nfeature {highly_variable_features} --resolution {params.paramaters[resolution]} "
        "--output.rds {output.rds} --output.xlsx {output.xlsx} --output.pca {output.pca} {umap_plot} {tsne_plot}"


rule umap_plot:
    input:
        "analyses/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds"
    output:
        umap="results/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/umap.plot.pdf"
    params:
        paramaters=paramspace.instance,
    shell:
        "workflow/scripts/scrna-umap.R --rds {input} --output.umap.plot {output.umap} --resolution {params.paramaters[resolution]}"


    
rule clustermarkers:
    input:
        "analyses/processed/{res}/{sample}.rds"
    output:
        "results/{sample}/resolution-{res}/{sample}.positive-markers-forAllClusters.xlsx",
        "results/{sample}/resolution-{res}/{sample}.all-markers-forAllClusters.xlsx"
    shell:
        """
        workflow/scripts/scrna-find-markers.R --rds {input} --resolution {wildcards.res} --sampleid {wildcards.sample} --logfc.threshold {logfc_threshold} --test.use {test_use}
        """

rule selected_markers:
    input:
        "analyses/processed/{res}/{sample}.rds"
    output:
        "results/{sample}/resolution-{res}/selected-markers/selected-markers-dotplot.pdf",
        directory("results/{sample}/resolution-{res}/selected-markers/plots/")
    shell:
        "workflow/scripts/scrna-selected-marker-plots.R --rds {input} --resolution {wildcards.res} --sampleid {wildcards.sample}"


rule positive_markers:
    input:
        rds="analyses/processed/{res}/{sample}.rds",
        excel="results/{sample}/resolution-{res}/{sample}.positive-markers-forAllClusters.xlsx"
    output:
        directory("results/{sample}/resolution-{res}/markers/")
    shell:
        "workflow/scripts/scrna-marker-plots.R --rds {input.rds} --resolution {wildcards.res} --sampleid {wildcards.sample} --xlsx {input.excel}"

rule integration_with_harmony:
    input:
        expand("analyses/processed/" + integration_resolution + "/{sample}.rds",sample=files)

    output:
        "analyses/integration/harmony/" + integration_id + "_harmony.rds"
    shell:
        """
        workflow/scripts/scrna-harmony.R --rds "{input}" --sampleid {integration_id} --resolution {integration_resolution} --normalization.method {normalization_method} --scale.factor {scale_factor} --nfeature {highly_variable_features}
        """


rule integration_with_seurat:
    input:
        expand("analyses/processed/" + integration_resolution + "/{sample}.rds",sample=files)
    output:
        "analyses/integration/seurat/" + integration_id + "_seurat.rds"
    shell:
        """
        workflow/scripts/scrna-seurat-integration.R --rds "{input}" --sampleid {integration_id} --resolution {integration_resolution}
        """


rule h5ad:
    input:
        rds="analyses/processed/{res}/{sample}.rds"
    output:
        "analyses/h5ad/{res}/{sample}.h5ad"
    shell:
        "workflow/scripts/scrna-convert-to-h5ad.R --rds {input.rds} --resolution {wildcards.res} --sampleid {wildcards.sample}"


rule celltype:
    input:
        "analyses/h5ad/{res}/{sample}.h5ad"
    
    output:
        directory("analyses/celltypist/{res}/{sample}/"),
        "analyses/celltypist/{res}/{sample}/predicted_labels.csv"

    shell:
        """
        celltypist --indata {input} --model {celltypist_model} --outdir {output[0]}
        """

rule seurat_celltype:
    input:
        csv="analyses/celltypist/{res}/{sample}/predicted_labels.csv",
        rds="analyses/processed/{res}/{sample}.rds"
    output:
        "analyses/celltypist/{res}/{sample}.rds"
    shell:
        """
        workflow/scripts/scrna-celltypist.R --rds {input.rds} --sampleid {wildcards.sample} --csv {input.csv} --output {output}
        """

rule go_enrichment:
    input:
        "results/{sample}/resolution-{res}/{sample}.all-markers-forAllClusters.xlsx"
    output:
        "results/{sample}/resolution-{res}/enrichment/GO-enrichment-all_clusters-ontology-{ontology}.xlsx"
    shell:
        """
        workflow/scripts/scrna-go_enrichment.R --xlsx {input} --output {output} --ontology {ontology} --algorithm {algorithm} --mapping {mapping} --statistics {statistics}
        """

