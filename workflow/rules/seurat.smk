from collections import defaultdict
from yaml import load
import os


def input_function(wildcards):
    if os.path.isfile(datafolder):
        return(datafolder)
    else:
        if os.path.isfile(datafolder + "/" + wildcards.sample + "/filtered_feature_bc_matrix/matrix.mtx.gz"):
            return(datafolder +  "/" + wildcards.sample + "/filtered_feature_bc_matrix/")
        elif os.path.isfile(datafolder + "/" + wildcards.sample + "/filtered_feature_bc_matrix.h5"):
            return(datafolder + "/" + wildcards.sample + "/")
        elif os.path.isfile(datafolder + "/" + wildcards.sample + "/raw_feature_bc_matrix/matrix.mtx.gz"):
            return(datafolder + "/" + wildcards.sample + "/raw_feature_bc_matrix/")
        elif os.path.isfile(datafolder + "/" + wildcards.sample + "/matrix.mtx.gz"):
            return(datafolder + "/" + wildcards.sample + "/")
        elif os.path.isfile(datafolder + "/" + wildcards.sample + "/matrix.mtx"):
            return(datafolder + "/" + wildcards.sample + "/")
        elif os.path.isfile(datafolder + "/" + wildcards.sample + "/outs/filtered_feature_bc_matrix/matrix.mtx.gz"):
            return(datafolder + "/" + wildcards.sample + "/outs/filtered_feature_bc_matrix/")
        else:
            return(datafolder + "/" + wildcards.sample + "/outs/raw_feature_bc_matrix/")

rule create_initial_raw_rds_and_trimming:
    input:
        #"data/{sample}/raw_feature_bc_matrix/"
        raw=input_function    
    output:
        rds=analyses_folder + "/raw/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds",
        before=[results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/technicals/before-qc-trimming-violinplot.pdf"] if is_integrated_sample is False else [],
        after=[results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/technicals/after-qc-trimming-violinplot.pdf"] if is_integrated_sample is False else [],
        mtplot=[results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/technicals/model-metrics-mitochondrial-genes.pdf"] if percent_mt == "auto" and is_integrated_sample is False else []
    params:
        mt_param=" --plot.mtplot " + results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/technicals/model-metrics-mitochondrial-genes.pdf" if percent_mt == "auto" and is_integrated_sample is False else " "
    

    run:
        if is_integrated_sample is True:
            shell("cp {input.raw} {output.rds}")
        else:
            shell("{cellsnake_path}workflow/scripts/scrna-read-qc.R --data.dir {input.raw} --output.rds {output.rds} --sampleid {wildcards.sample} --percent.rp {percent_rp} --percent.mt {wildcards.percent_mt} --min.features {min_features} --min.cells {min_cells} --before.violin.plot {output.before} --after.violin.plot {output.after} {params.mt_param}")
            


rule clustree:
    input:
        analyses_folder + "/raw/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds"
    output:
        clustree=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/ClustTree.pdf",
        heatmap=[results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/technicals/DimHeatMap_plot.pdf"] if is_integrated_sample is False else [],
        hvfplot=[results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/technicals/highly-variable-features.pdf"] if is_integrated_sample is False else [],
        jackandelbow=[results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/technicals/JackandElbow_plot.pdf"] if is_integrated_sample is False else []
    run:
        if is_integrated_sample is True:
            shell("{cellsnake_path}workflow/scripts/scrna-clusteringtree.R --rds {input} --output {output.clustree} --integration")
        else:
            shell("{cellsnake_path}workflow/scripts/scrna-clusteringtree.R --rds {input} --output {output.clustree} --heatmap {output.heatmap} --hvfplot {output.hvfplot} --jackandelbow {output.jackandelbow}")



rule normalization_pca_rds:
    input:
        analyses_folder + "/raw/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds"
    output:
        rds=analyses_folder + "/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds",
        xlsx=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/number-of-cells-per-cluster.xlsx"
    params:
        paramaters=paramspace.instance,
        doublet_filter=doublet_filter,
        umap_plot=umap_plot,
        tsne_plot=tsne_plot,
        integration="--integration" if is_integrated_sample is True else " "
    shell:
        "{cellsnake_path}workflow/scripts/scrna-normalization-pca.R --rds {input} {params.doublet_filter} --normalization.method {normalization_method} "
        "--scale.factor {scale_factor} --nfeature {highly_variable_features} --resolution {params.paramaters[resolution]} "
        "--output.rds {output.rds} --output.xlsx {output.xlsx} {umap_plot} {tsne_plot} {params.integration}"

rule dim_plots:
    input:
        analyses_folder + "/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds"
    output:
        pl=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/{reduction}-{i}.plot.pdf"
    params:
        paramaters=paramspace.instance,
    shell:
        "{cellsnake_path}workflow/scripts/scrna-dimplot.R --rds {input} --reduction.type {wildcards.reduction} --output.reduction.plot {output.pl} --idents {wildcards.i}"

    
    
rule find_all_cluster_markers:
    input:
        analyses_folder + "/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds"
    output:
        positive=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/positive-markers-forAllClusters.xlsx",
        allmarkers=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  +  "/all-markers-forAllClusters.xlsx"
    params:
        paramaters=paramspace.instance,
    shell:
        "{cellsnake_path}workflow/scripts/scrna-find-markers.R --rds {input} --resolution {params.paramaters[resolution]} --logfc.threshold {logfc_threshold} --test.use {test_use} --output.xlsx.positive {output.positive} --output.xlsx.all {output.allmarkers}"


rule plot_top_markers:
    input:
        results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/positive-markers-forAllClusters.xlsx"
    output:
        results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/summarized-cluster-markers-plot.pdf"
    shell:
        "{cellsnake_path}workflow/scripts/scrna-top-marker-plot.R --xlsx {input} --output.plot {output}" 


rule plot_top_positive_markers:
    input:
        rds=analyses_folder + "/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds",
        excel=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/positive-markers-forAllClusters.xlsx"
    output:
        directory(results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  +  "/positive_marker_plots_{reduction}/")
    params:
        paramaters=paramspace.instance,
    shell:
        "{cellsnake_path}workflow/scripts/scrna-marker-plots.R --rds {input.rds} --xlsx {input.excel} --top_n {marker_plots_per_cluster_n} --output.plot.dir {output} --reduction.type {wildcards.reduction}"



rule selected_marker_dot_plot:
    input:
        analyses_folder + "/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds"
    output:
        dotplot=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/selected-genes-dotplot.pdf"
    params:
        paramaters=paramspace.instance,
    shell:
        "{cellsnake_path}workflow/scripts/scrna-dotplot.R --rds {input} --tsv {selected_markers_file} --resolution {params.paramaters[resolution]} --output.dotplot {output.dotplot}"


rule selected_marker_plots_file:
    input:
        analyses_folder + "/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds"
    output:
        sdir=directory(results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  +  "/selected_gene_plots_from_file_{reduction}/")
    params:
        paramaters=paramspace.instance,
    shell:
        "{cellsnake_path}workflow/scripts/scrna-selected-marker-plots.R --rds {input} --tsv {selected_markers_file} --output.plot.dir {output.sdir} --reduction.type {wildcards.reduction}"


rule selected_marker_plot:
    input:
        analyses_folder + "/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds"
    output:
        out=[results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  +  "/selected_gene_plots_{reduction}/" + g + ".pdf" for g in gene_to_plot if g is not None]
    params:
        paramaters=paramspace.instance,
        sdir=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  +  "/selected_gene_plots_{reduction}/",
        gene=str(" ".join(gene_to_plot))
    shell:
        """{cellsnake_path}workflow/scripts/scrna-selected-marker-plots.R --rds {input} --gene "{params.gene}" --output.plot.dir {params.sdir} --reduction.type {wildcards.reduction}"""




rule h5ad:
    input:
        rds=analyses_folder + "/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds"
    output:
        analyses_folder + "/h5ad/" + f"{paramspace.wildcard_pattern}" + "/{sample}.h5ad"
    shell:
        "{cellsnake_path}workflow/scripts/scrna-convert-to-h5ad.R --rds {input.rds} --output {output}"


rule celltype:
    input:
        analyses_folder + "/h5ad/" + f"{paramspace.wildcard_pattern}" + "/{sample}.h5ad"
    
    output:
        outputdir=directory(analyses_folder + "/celltypist/" + celltypist_model + "/" + f"{paramspace.wildcard_pattern}" + "/{sample}"),
        predicted=analyses_folder + "/celltypist/" + celltypist_model + "/" + f"{paramspace.wildcard_pattern}" + "/{sample}/predicted_labels.csv",
        dotplot=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/celltype_annotation/" + celltypist_model + "/annotation.dotplot.pdf",
        xlsx=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/celltype_annotation/" + celltypist_model + "/cluster_annotation_table.xlsx"
        
    shell:
        """
        #celltypist --indata {input} --model {celltypist_model} --majority-voting --outdir {output[0]} 
        {cellsnake_path}workflow/scripts/scrna-celltypist.py {input} {output.dotplot} {output.outputdir} {output.xlsx} {celltypist_model}
        """

rule seurat_celltype:
    input:
        csv=analyses_folder + "/celltypist/" + celltypist_model + "/" + f"{paramspace.wildcard_pattern}" + "/{sample}/predicted_labels.csv",
        rds=analyses_folder + "/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds"
    output:
        umap=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/celltype_annotation/" + celltypist_model + "/annotation.umap.pdf",
        tsne=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/celltype_annotation/" + celltypist_model + "/annotation.tsne.pdf"
    shell:
        """
        {cellsnake_path}workflow/scripts/scrna-celltypist.R --rds {input.rds} --csv {input.csv} --output.tsne.plot {output.tsne} --output.umap.plot {output.umap}
        """

rule go_enrichment:
    input:
        results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  +  "/all-markers-forAllClusters.xlsx"
    output:
        results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/enrichment_analysis/GO-enrichment-" + ontology +  "-all_clusters.xlsx"
    shell:
        """
        {cellsnake_path}workflow/scripts/scrna-go_enrichment.R --xlsx {input} --output {output} --ontology {ontology} --algorithm {algorithm} --mapping {mapping} --statistics {statistics}
        """

rule kegg_enrichment:
    input:
        results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  +  "/all-markers-forAllClusters.xlsx"
    output:
        rds=analyses_folder + "/kegg/" + f"{paramspace.wildcard_pattern}" + "/{sample}_kegg.rds",
        kegg=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/enrichment_analysis/KEGG-over_representation-all_clusters.xlsx",
        gse=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/enrichment_analysis/KEGG-geneset-all_clusters.xlsx",
        mkegg=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/enrichment_analysis/KEGG-module_over_representation-all_clusters.xlsx",
        mgse=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/enrichment_analysis/KEGG-module_geneset-all_clusters.xlsx"
    shell:
        "{cellsnake_path}workflow/scripts/scrna-kegg.R --xlsx {input} --output.rds {output.rds} --mapping {mapping} --output.kegg {output.kegg} --output.mkegg {output.mkegg}  --output.gse {output.gse} --output.mgse {output.mgse}"



rule gsea:
    input:
        rds=analyses_folder + "/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds",
        csv=lambda wildcards: analyses_folder + "/celltypist/" + celltypist_model + "/" + f"{paramspace.wildcard_pattern}" + "/{sample}/predicted_labels.csv" if wildcards.i == "majority_voting" else [],
        gseafile=gsea_file
    output:
        xlsx=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/gsea/gsea-{i}-output.xlsx"

    run:
        if wildcards.i == "majority_voting":
            shell("{cellsnake_path}workflow/scripts/scrna-gsea.R --idents {wildcards.i} --rds {input.rds} --gseafile {input.gseafile} --output.xlsx {output.xlsx}  --csv {input.csv}")
        else:
            shell("{cellsnake_path}workflow/scripts/scrna-gsea.R --idents {wildcards.i} --rds {input.rds} --gseafile {input.gseafile} --output.xlsx {output.xlsx}")


rule cellchat:
    input:
        rds=analyses_folder + "/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds",
        csv=lambda wildcards: analyses_folder + "/celltypist/" + celltypist_model + "/" + f"{paramspace.wildcard_pattern}" + "/{sample}/predicted_labels.csv" if wildcards.i == "majority_voting" else []
    output:
        cellchatrds=analyses_folder + "/cellchat/" + f"{paramspace.wildcard_pattern}" + "/{sample}/cellchat_{i}.rds"
    threads: 5
    run:
        if wildcards.i == "majority_voting":
            shell("{cellsnake_path}workflow/scripts/scrna-cellchat.R --rds {input.rds} --species {species} --idents {wildcards.i} --output {output.cellchatrds} --csv {input.csv}")
        else:
            shell("{cellsnake_path}workflow/scripts/scrna-cellchat.R --rds {input.rds} --species {species} --idents {wildcards.i} --output {output.cellchatrds}")

rule cellchat_plots:
    input:
        cellchatrds=analyses_folder + "/cellchat/" + f"{paramspace.wildcard_pattern}" + "/{sample}/cellchat_{i}.rds"
    output:
        outputdir=directory(results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/cellchat/{i}/")
    shell:
        "{cellsnake_path}workflow/scripts/scrna-cellchat_plots.R --rds {input.cellchatrds} --output.dir {output.outputdir}"
