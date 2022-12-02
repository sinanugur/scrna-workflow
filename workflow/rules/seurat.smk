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


rule create_initial_raw_rds_and_trimming:
    input:
        #"data/{sample}/raw_feature_bc_matrix/"
        raw=input_function    
    output:
        rds="analyses/raw/{sample}/" + f"{paramspace.wildcard_pattern}.rds",
        before=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/technicals/before-qc-trimming-violinplot.pdf",
        after=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/technicals/after-qc-trimming-violinplot.pdf",
        mtplot=[results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/technicals/model-metrics-mitochondrial-genes.pdf"] if auto_mt_filtering else []
    run:
        if auto_mt_filtering:
            shell("{cellsnake_path}workflow/scripts/scrna-read-qc.R --data.dir {input.raw} --output.rds {output.rds} --sampleid {wildcards.sample} --percent.rp {percent_rp} --min.features {min_features} --min.cells {min_cells} --auto.mt.filter --plot.mtplot {output.mtplot} --before.violin.plot {output.before} --after.violin.plot {output.after}")
        else:
            shell("{cellsnake_path}workflow/scripts/scrna-read-qc.R --data.dir {input.raw} --output.rds {output.rds} --sampleid {wildcards.sample} --percent.rp {percent_rp} --percent.mt {wildcards.MT} --min.features {min_features} --min.cells {min_cells} --before.violin.plot {output.before} --after.violin.plot {output.after}")
            


rule clustree:
    input:
        "analyses/raw/{sample}/" + f"{paramspace.wildcard_pattern}.rds"
    output:
        clustree=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/clusteringTree/clusteringTree.pdf",
        heatmap=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/technicals/DimHeatMap_plot.pdf",
        hvfplot=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/technicals/highly-variable-features.pdf",
        jackandelbow=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/technicals/JackandElbow_plot.pdf"
    shell:
        "{cellsnake_path}workflow/scripts/scrna-clusteringtree.R --rds {input} --output {output.clustree} --heatmap {output.heatmap} --hvfplot {output.hvfplot} --jackandelbow {output.jackandelbow}"



rule normalization_pca_rds:
    input:
        "analyses/raw/{sample}/" + f"{paramspace.wildcard_pattern}.rds"
    output:
        rds="analyses/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds",
        xlsx=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/number-of-cells-per-cluster.xlsx",
        pca=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/pca.plot.pdf"
    params:
        paramaters=paramspace.instance,
        doublet_filter=doublet_filter,
        umap_plot=umap_plot,
        tsne_plot=tsne_plot
    shell:
        "{cellsnake_path}workflow/scripts/scrna-normalization-pca.R --rds {input} {params.doublet_filter} --normalization.method {normalization_method} "
        "--scale.factor {scale_factor} --nfeature {highly_variable_features} --resolution {params.paramaters[resolution]} "
        "--output.rds {output.rds} --output.xlsx {output.xlsx} --output.pca {output.pca} {umap_plot} {tsne_plot}"


rule umap_plot:
    input:
        "analyses/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds"
    output:
        umap=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/umap.plot.pdf"
    params:
        paramaters=paramspace.instance,
    shell:
        "{cellsnake_path}workflow/scripts/scrna-dimplot.R --rds {input} --reduction.type umap --output.reduction.plot {output.umap} --resolution {params.paramaters[resolution]}"

rule tsne_plot:
    input:
        "analyses/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds"
    output:
        tsne=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/tsne.plot.pdf"
    params:
        paramaters=paramspace.instance,
    shell:
        "{cellsnake_path}workflow/scripts/scrna-dimplot.R --rds {input} --reduction.type tsne --output.reduction.plot {output.tsne} --resolution {params.paramaters[resolution]}"


    
rule find_all_cluster_markers:
    input:
        "analyses/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds"
    output:
        positive=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/positive-markers-forAllClusters.xlsx",
        allmarkers=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  +  "/all-markers-forAllClusters.xlsx"
    params:
        paramaters=paramspace.instance,
    shell:
        "{cellsnake_path}workflow/scripts/scrna-find-markers.R --rds {input} --resolution {params.paramaters[resolution]} --logfc.threshold {logfc_threshold} --test.use {test_use} --output.xlsx.positive {output.positive} --output.xlsx.all {output.allmarkers}"


rule plot_top_positive_markers:
    input:
        rds="analyses/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds",
        excel=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/positive-markers-forAllClusters.xlsx"
    output:
        directory(results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  +  "/positive_marker_plots_{reduction}/")
    params:
        paramaters=paramspace.instance,
    shell:
        "{cellsnake_path}workflow/scripts/scrna-marker-plots.R --rds {input.rds} --resolution {params.paramaters[resolution]} --xlsx {input.excel} --top_n {marker_plots_per_cluster_n} --output.plot.dir {output} --reduction.type {wildcards.reduction}"



rule selected_marker_dot_plot:
    input:
        "analyses/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds"
    output:
        dotplot=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/selected-markers-dotplot.pdf"
    params:
        paramaters=paramspace.instance,
    shell:
        "{cellsnake_path}workflow/scripts/scrna-dotplot.R --rds {input} --tsv {selected_markers_file} --resolution {params.paramaters[resolution]} --output.dotplot {output.dotplot}"


rule selected_marker_plots:
    input:
        "analyses/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds"
    output:
        sdir=directory(results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  +  "/selected_marker_plots_{reduction}/")
    params:
        paramaters=paramspace.instance,
    shell:
        "{cellsnake_path}workflow/scripts/scrna-selected-marker-plots.R --rds {input} --tsv {selected_markers_file} --resolution {params.paramaters[resolution]} --output.plot.dir {output.sdir} --reduction.type {wildcards.reduction}"





rule integration_with_harmony:
    input:
        expand("analyses/processed/" + integration_resolution + "/{sample}.rds",sample=files)

    output:
        "analyses/integration/harmony/" + integration_id + "_harmony.rds"
    shell:
        """
        {cellsnake_path}workflow/scripts/scrna-harmony.R --rds "{input}" --sampleid {integration_id} --resolution {integration_resolution} --normalization.method {normalization_method} --scale.factor {scale_factor} --nfeature {highly_variable_features}
        """


rule integration_with_seurat:
    input:
        expand("analyses/processed/" + integration_resolution + "/{sample}.rds",sample=files)
    output:
        "analyses/integration/seurat/" + integration_id + "_seurat.rds"
    shell:
        """
        {cellsnake_path}workflow/scripts/scrna-seurat-integration.R --rds "{input}" --sampleid {integration_id} --resolution {integration_resolution}
        """


rule h5ad:
    input:
        rds="analyses/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds"
    output:
        "analyses/h5ad/" + f"{paramspace.wildcard_pattern}" + "/{sample}.h5ad"
    shell:
        "{cellsnake_path}workflow/scripts/scrna-convert-to-h5ad.R --rds {input.rds} --output {output}"


rule celltype:
    input:
        "analyses/h5ad/" + f"{paramspace.wildcard_pattern}" + "/{sample}.h5ad"
    
    output:
        outputdir=directory("analyses/celltypist/" + f"{paramspace.wildcard_pattern}" + "/{sample}"),
        predicted="analyses/celltypist/" + f"{paramspace.wildcard_pattern}" + "/{sample}/predicted_labels.csv",
        dotplot=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/celltype_annotation/annotation.dotplot.pdf",
        xlsx=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/celltype_annotation/cluster_annotation_table.xlsx"
        
    shell:
        """
        #celltypist --indata {input} --model {celltypist_model} --majority-voting --outdir {output[0]} 
        {cellsnake_path}workflow/scripts/scrna-celltypist.py {input} {output.dotplot} {output.outputdir} {output.xlsx} {celltypist_model}
        """

rule seurat_celltype:
    input:
        csv="analyses/celltypist/" + f"{paramspace.wildcard_pattern}" + "/{sample}/predicted_labels.csv",
        rds="analyses/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds"
    output:
        umap=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/celltype_annotation/annotation.umap.pdf",
        tsne=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/celltype_annotation/annotation.tsne.pdf"
    shell:
        """
        {cellsnake_path}workflow/scripts/scrna-celltypist.R --rds {input.rds} --csv {input.csv} --output.tsne.plot {output.tsne} --output.umap.plot {output.umap}
        """

rule go_enrichment:
    input:
        results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  +  "/all-markers-forAllClusters.xlsx"
    output:
        results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "GO-enrichment-" + ontology +  "-all_clusters.xlsx"
    shell:
        """
        {cellsnake_path}workflow/scripts/scrna-go_enrichment.R --xlsx {input} --output {output} --ontology {ontology} --algorithm {algorithm} --mapping {mapping} --statistics {statistics}
        """

rule gsea:
    input:
        rds="analyses/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds",
        gseafile=gsea_file
    output:
        directory(results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/gsea/")
    shell:
        """
        {cellsnake_path}workflow/scripts/scrna-gsea.R --rds {input.rds} --gseafile {input.gseafile} --output.dir {output}
        """