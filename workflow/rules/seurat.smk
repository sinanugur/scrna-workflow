from collections import defaultdict
from yaml import load
import os


def input_function(wildcards):
    if os.path.isfile("data/" + wildcards.sample + "/raw_feature_bc_matrix/matrix.mtx.gz"):
        return("data/" + wildcards.sample + "/raw_feature_bc_matrix/")
    else:
        return("data/" + wildcards.sample + "/outs/filtered_feature_bc_matrix/")

rule rds:
    input:
        #"data/{sample}/raw_feature_bc_matrix/"
        input_function
    output:
        "analyses/rawrds/{sample}.rds"

    shell:
        "workflow/scripts/scrna-read-qc.R --data.dir {input} --sampleid {wildcards.sample} --percent.mt {percent_mt} --min.features {min_features} --min.cells {min_cells} --minCov {min_coverage}"

rule process_rds:
    input:
        "analyses/rawrds/{sample}.rds"
    output:
        "analyses/processed/{sample}.rds"
    shell:
        "workflow/scripts/scrna-normalization-pca.R --rds {input} --sampleid {wildcards.sample} --normalization.method {normalization_method} --scale.factor {scale_factor} --nfeature {highly_variable_features}"


rule clustree:
    input:
        "analyses/processed/{sample}.rds"
    output:
        "results/{sample}/clusteringTree/clusteringTree-{sample}.pdf"
    shell:
        "workflow/scripts/scrna-clusteringtree.R --rds {input} --sampleid {wildcards.sample}"
    
rule clustermarkers:
    input:
        "analyses/processed/{sample}.rds"
    output:
        "results/{sample}/resolution-{res}/{sample}.positive-markers-forAllClusters.xlsx",
        "results/{sample}/resolution-{res}/{sample}.all-markers-forAllClusters.xlsx"
    shell:
        """
        workflow/scripts/scrna-find-markers.R --rds {input} --resolution {wildcards.res} --sampleid {wildcards.sample} --logfc.threshold {logfc_threshold} --test.use {test_use}
        """

rule metrics:
    input:
        "analyses/processed/{sample}.rds"
    output:
        "results/{sample}/resolution-{res}/{sample}.number-of-cells-per-cluster.xlsx"
    shell:
        "workflow/scripts/scrna-metrics.R --rds {input} --resolution {wildcards.res} --sampleid {wildcards.sample}"

rule selected_markers:
    input:
        "analyses/processed/{sample}.rds"
    output:
        "results/{sample}/resolution-{res}/selected-markers-dotplot.pdf",
        directory("results/{sample}/resolution-{res}/selected-markers/plots/")
    shell:
        "workflow/scripts/scrna-selected-marker-plots.R --rds {input} --resolution {wildcards.res} --sampleid {wildcards.sample}"


rule positive_markers:
    input:
        rds="analyses/processed/{sample}.rds",
        excel="results/{sample}/resolution-{res}/{sample}.positive-markers-forAllClusters.xlsx"
    output:
        directory("results/{sample}/resolution-{res}/markers/")
    shell:
        "workflow/scripts/scrna-marker-plots.R --rds {input.rds} --resolution {wildcards.res} --sampleid {wildcards.sample} --xlsx {input.excel}"

rule integration_with_harmony:
    input:
        expand("analyses/processed/{sample}.rds",sample=files)
    output:
        "analyses/integration/harmony/" + integration_id + "_harmony.rds"
    shell:
        """
        workflow/scripts/scrna-harmony.R --rds "{input}" --sampleid {integration_id} --resolution {integration_resolution} --normalization.method {normalization_method} --scale.factor {scale_factor} --nfeature {highly_variable_features}
        """


rule integration_with_seurat:
    input:
        expand("analyses/processed/{sample}.rds",sample=files)
    output:
        "analyses/integration/seurat/" + integration_id + "_seurat.rds"
    shell:
        """
        workflow/scripts/scrna-seurat-integration.R --rds "{input}" --sampleid {integration_id} --resolution {integration_resolution}
        """
