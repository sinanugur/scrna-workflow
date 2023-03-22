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
        before=[results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/technicals/plot_before-qc-trimming.pdf"] if is_integrated_sample is False else [],
        after=[results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/technicals/plot_after-qc-trimming.pdf"] if is_integrated_sample is False else [],
        mtplot=[results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/technicals/plot_model-metrics-mitochondrial-genes.pdf"] if percent_mt == "auto" and is_integrated_sample is False else []
    params:
        mt_param=" --plot.mtplot " + results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/technicals/plot_model-metrics-mitochondrial-genes.pdf" if percent_mt == "auto" and is_integrated_sample is False else " "
    

    run:
        if is_integrated_sample is True:
            shell("cp {input.raw} {output.rds}")
        else:
            shell("{cellsnake_path}workflow/scripts/scrna-read-qc.R --data.dir {input.raw} --output.rds {output.rds} --sampleid {wildcards.sample} --percent.rp {percent_rp} --percent.mt {wildcards.percent_mt} --min.features {min_features} --max.features {max_features} --max.molecules {max_molecules}  --min.cells {min_cells} --before.violin.plot {output.before} --after.violin.plot {output.after} {params.mt_param}")
            


rule clustree:
    input:
        analyses_folder + "/raw/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds"
    output:
        clustree=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/technicals/plot_clustree.pdf",
        heatmap=[results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/technicals/plot_dimheatmap.pdf"] if is_integrated_sample is False else [],
        hvfplot=[results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/technicals/plot_highly_variable_features.pdf"] if is_integrated_sample is False else [],
        jackandelbow=[results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/technicals/plot_jackandelbow.pdf"] if is_integrated_sample is False else []
    run:
        if is_integrated_sample is True:
            shell("{cellsnake_path}workflow/scripts/scrna-clusteringtree.R --rds {input} --clplot {output.clustree} --integration")
        else:
            shell("{cellsnake_path}workflow/scripts/scrna-clusteringtree.R --rds {input} --scale.factor {scale_factor} --nfeature {highly_variable_features} --variable.selection.method {variable_selection_method} --normalization.method {normalization_method} --clplot {output.clustree} --heplot {output.heatmap} --hvfplot {output.hvfplot} --jeplot {output.jackandelbow}")



rule normalization_pca_rds:
    input:
        analyses_folder + "/raw/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds"
    output:
        rds=analyses_folder + "/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds"
    params:
        paramaters=paramspace.instance,
        doublet_filter=doublet_filter,
        umap_plot=umap_plot,
        tsne_plot=tsne_plot,
        integration="--integration" if is_integrated_sample is True else " "
    shell:
        "{cellsnake_path}workflow/scripts/scrna-normalization-pca.R --rds {input} {params.doublet_filter} --normalization.method {normalization_method} "
        "--scale.factor {scale_factor} --reference {singler_ref} --variable.selection.method {variable_selection_method} --nfeature {highly_variable_features} --resolution {params.paramaters[resolution]} "
        "--output.rds {output.rds} {umap_plot} {tsne_plot} {params.integration}"

rule some_metrics:
    input:
        rds=analyses_folder + "/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds"
    output:
        ccplot=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/metrics/plot_cellcount-{i}.pdf",
        ccbarplot=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/metrics/plot_cellcount_barplot-{i}.pdf",
        html=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/metrics/plot_cellcount_barplot-{i}.html",
        t=temp(directory(results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/metrics/plot_cellcount_barplot-{i}_files/"))
    shell:
        "{cellsnake_path}workflow/scripts/scrna-metrics.R --rds {input.rds} --ccplot {output.ccplot} --ccbarplot {output.ccbarplot} --html {output.html} --idents {wildcards.i}"

rule some_technicals:
    input:
        rds=analyses_folder + "/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds"
    output:
        fplot=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/technicals/plot_nFeature.pdf",
        cplot=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/technicals/plot_nCount.pdf"

    shell:
        "{cellsnake_path}workflow/scripts/scrna-technicals.R --rds {input.rds} --fplot {output.fplot} --cplot {output.cplot}"



rule dim_plots:
    input:
        analyses_folder + "/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds"
    output:
        pl=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/plot_dimplot_{reduction}-{i}.pdf",
        html=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/plot_dimplot_{reduction}-{i}.html",
        f=temp(directory(results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/plot_dimplot_{reduction}-{i}_files/"))
    params:
        paramaters=paramspace.instance,
    shell:
        "{cellsnake_path}workflow/scripts/scrna-dimplot.R --rds {input} --reduction.type {wildcards.reduction} --pdfplot {output.pl} --htmlplot {output.html} --idents {wildcards.i}"



rule singler_dim_plots: #this is singler annotation that does not include cluster idents and only use singler annotation
    input:
        analyses_folder + "/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds"
    output:
        pl=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/plot_annotation_{reduction}.pdf",
        html=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/plot_annotation_{reduction}.html",
        t=temp(directory(results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/plot_annotation_{reduction}_files"))
    params:
        paramaters=paramspace.instance,
    shell:
        "{cellsnake_path}workflow/scripts/scrna-dimplot.R --rds {input} --reduction.type {wildcards.reduction} --pdfplot {output.pl} --htmlplot {output.html} --idents singler"

    
    
rule find_all_cluster_markers:
    input:
        analyses_folder + "/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds"
    output:
        positive=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/table_positive-markers-{i}.xlsx",
        allmarkers=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  +  "/table_all-markers-{i}.xlsx"
    params:
        paramaters=paramspace.instance,
    shell:
        "{cellsnake_path}workflow/scripts/scrna-find-markers.R --rds {input} --idents {wildcards.i} --resolution {params.paramaters[resolution]} --logfc.threshold {logfc_threshold} --test.use {test_use} --output.xlsx.positive {output.positive} --output.xlsx.all {output.allmarkers}"


rule plot_top_markers:
    input:
        results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/table_positive-markers-{i}.xlsx"
    output:
        results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/summarized_markers-for-{i}.pdf"
    shell:
        "{cellsnake_path}workflow/scripts/scrna-top-marker-plot.R --xlsx {input} --output.plot {output}" 


rule plot_top_positive_markers:
    input:
        rds=analyses_folder + "/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds",
        excel=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  + "/table_positive-markers-{i}.xlsx"
    output:
        directory(results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  +  "/positive_marker_plots_{reduction}/{i}/")
    params:
        paramaters=paramspace.instance,
    shell:
        "{cellsnake_path}workflow/scripts/scrna-marker-plots.R --rds {input.rds} --xlsx {input.excel} --top_n {marker_plots_per_cluster_n} --output.plot.dir {output} --reduction.type {wildcards.reduction}"



rule selected_marker_plot:
    input:
        analyses_folder + "/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds"
    output:
        out=[results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  +  "/selected_gene_plots_{reduction}/{i}/" + g + ".pdf" for g in gene_to_plot if g is not None]
    params:
        paramaters=paramspace.instance,
        sdir=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  +  "/selected_gene_plots_{reduction}/{i}/",
        gene=str(" ".join(gene_to_plot))
    shell:
        """{cellsnake_path}workflow/scripts/scrna-selected-marker-plots.R --rds {input} --gene "{params.gene}" --output.plot.dir {params.sdir} --reduction.type {wildcards.reduction} --idents {wildcards.i}"""




rule h5ad:
    input:
        rds=analyses_folder + "/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds"
    output:
        analyses_folder + "/h5ad/" + f"{paramspace.wildcard_pattern}" + "/{sample}.h5ad"
    shell:
        "{cellsnake_path}workflow/scripts/scrna-convert-to-h5ad.R --rds {input.rds} --output {output}"


rule singler_celltype: #this singler rule use idents information from the clustering step
    input:
        rds=analyses_folder + "/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds"
    output:
        sheplot=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/singler/plot_score_heatmap-{i}.pdf",
        pheplot=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/singler/plot_clusters-{i}.pdf",
        sheplottop=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/singler/plot_score_heatmap_top-{i}.pdf"

    shell:
        "{cellsnake_path}workflow/scripts/scrna-singler-plots.R --rds {input.rds} --reference {singler_ref} --sheplot {output.sheplot} --pheplot {output.pheplot} --sheplottop {output.sheplottop} --idents {wildcards.i}"


rule celltypist_celltype:
    input:
        analyses_folder + "/h5ad/" + f"{paramspace.wildcard_pattern}" + "/{sample}.h5ad"
    
    output:
        outputdir=directory(analyses_folder + "/celltypist/" + celltypist_model + "/" + f"{paramspace.wildcard_pattern}" + "/{sample}/{i}/"),
        predicted=analyses_folder + "/celltypist/" + celltypist_model + "/" + f"{paramspace.wildcard_pattern}" + "/{sample}/{i}/predicted_labels.csv",
        dotplot=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/celltypist/" + celltypist_model + "/plot_celltypist_dotplot-{i}.pdf",
        xlsx=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/celltypist/" + celltypist_model + "/table_cluster_annotation_table-{i}.xlsx"
        
    shell:
        "{cellsnake_path}workflow/scripts/scrna-celltypist.py {input} {output.dotplot} {output.outputdir} {output.xlsx} {celltypist_model} {wildcards.i}"

rule seurat_celltype_celltypist:
    input:
        csv=analyses_folder + "/celltypist/" + celltypist_model + "/" + f"{paramspace.wildcard_pattern}" + "/{sample}/{i}/predicted_labels.csv",
        rds=analyses_folder + "/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds"
    output:
        umap=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/celltypist/" + celltypist_model + "/plot_celltypist_umap-{i}.pdf",
        tsne=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/celltypist/" + celltypist_model + "/plot_celltypist_tsne-{i}.pdf"
    shell:
        """
        {cellsnake_path}workflow/scripts/scrna-celltypist.R --rds {input.rds} --csv {input.csv} --output.tsne.plot {output.tsne} --output.umap.plot {output.umap}
        """


rule kegg_enrichment:
    input:
        results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  +  "/table_all-markers-{i}.xlsx"
    output:
        rds=analyses_folder + "/kegg/" + f"{paramspace.wildcard_pattern}" + "/{sample}-{i}-kegg.rds",
        kegg=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/enrichment_analysis/table_KEGG-enrichment-{i}.xlsx",
        gse=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/enrichment_analysis/table_KEGG-geneset_enrichment-{i}.xlsx",
        mkegg=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/enrichment_analysis/table_KEGG-module_enrichment-{i}.xlsx",
        mgse=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/enrichment_analysis/table_KEGG-module_geneset_enrichment-{i}.xlsx"
    shell:
        "{cellsnake_path}workflow/scripts/scrna-kegg.R --xlsx {input} --output.rds {output.rds} --mapping {mapping} --organism {organism} --output.kegg {output.kegg} --output.mkegg {output.mkegg}  --output.gse {output.gse} --output.mgse {output.mgse}"

rule go2_enrichment:
    input:
        results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}"  +  "/table_all-markers-{i}.xlsx"
    output:
        rds=analyses_folder + "/go/" + f"{paramspace.wildcard_pattern}" + "/{sample}-{i}-go.rds",
        go=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/enrichment_analysis/table_GO-enrichment-{i}.xlsx",
        gse=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/enrichment_analysis/table_GO-geneset_enrichment-{i}.xlsx"
    shell:
        "{cellsnake_path}workflow/scripts/scrna-go_analysis.R --xlsx {input} --output.rds {output.rds} --mapping {mapping} --output.go {output.go} --output.gse {output.gse}"



rule gsea_cerebro:
    input:
        rds=analyses_folder + "/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds",
        #csv=lambda wildcards: analyses_folder + "/celltypist/" + celltypist_model + "/" + f"{paramspace.wildcard_pattern}" + "/{sample}/predicted_labels.csv" if wildcards.i == "majority_voting" else [],
        gseafile=gsea_file
    output:
        xlsx=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/gsea/table_gsea_output-{i}.xlsx"

    run:
        if wildcards.i == "majority_voting":
            shell("{cellsnake_path}workflow/scripts/scrna-gsea.R --idents {wildcards.i} --rds {input.rds} --gseafile {input.gseafile} --output.xlsx {output.xlsx}  --csv {input.csv}")
        else:
            shell("{cellsnake_path}workflow/scripts/scrna-gsea.R --idents {wildcards.i} --rds {input.rds} --gseafile {input.gseafile} --output.xlsx {output.xlsx}")


rule cellchat:
    input:
        rds=analyses_folder + "/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds"
        #csv=lambda wildcards: analyses_folder + "/celltypist/" + celltypist_model + "/" + f"{paramspace.wildcard_pattern}" + "/{sample}/predicted_labels.csv" if wildcards.i == "majority_voting" else []
    output:
        cellchatrds=analyses_folder + "/cellchat/" + f"{paramspace.wildcard_pattern}" + "/{sample}/cellchat-{i}.rds"
    threads: 5
    run:
        if wildcards.i == "majority_voting":
            shell("{cellsnake_path}workflow/scripts/scrna-cellchat.R --rds {input.rds} --species {species} --idents {wildcards.i} --output {output.cellchatrds} --csv {input.csv}")
        else:
            shell("{cellsnake_path}workflow/scripts/scrna-cellchat.R --rds {input.rds} --species {species} --idents {wildcards.i} --output {output.cellchatrds}")

rule cellchat_plots:
    input:
        cellchatrds=analyses_folder + "/cellchat/" + f"{paramspace.wildcard_pattern}" + "/{sample}/cellchat-{i}.rds"
    output:
        outputdir=directory(results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/cellchat/{i}/")
    shell:
        "{cellsnake_path}workflow/scripts/scrna-cellchat_plots.R --rds {input.cellchatrds} --output.dir {output.outputdir}"

rule monocle3_plots:
    input:
        rds=analyses_folder + "/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds"
    output:
        pplot=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/trajectory/plot_monocle-partition-plot.pdf"
    params:
        outputdir=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/trajectory/"
    shell:
        "{cellsnake_path}workflow/scripts/scrna-monocle3.R --rds {input.rds} --output.dir {params.outputdir} --pplot {output.pplot}"

rule subset_final_rds:
    input:
        rds=analyses_folder + "/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds"
    output:
        output=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/dataoutput/{keyword}.rds"
    shell:
        "{cellsnake_path}workflow/scripts/scrna-subset_final_rds.R --rds {input.rds} --output {output.rds}"