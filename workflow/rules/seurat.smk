from collections import defaultdict
from yaml import load


rule rds:
    input:
        "data/{sample}/raw_feature_bc_matrix/"
    output:
        "analyses/{sample}/rawrds/{sample}.rds"

    shell:
        "workflow/scripts/scrna-read-qc.R --data.dir {input} --sampleid {wildcards.sample} --percent.mt {percent_mt} --min.features {min_features} --min.cells {min_cells} --minCov {min_coverage}"

rule process_rds:
    input:
        "analyses/{sample}/rawrds/{sample}.rds"
    output:
        "analyses/{sample}/processed/{sample}.rds"
    shell:
        "workflow/scripts/scrna-normalization-pca.R --rds {input} --sampleid {wildcards.sample}"


rule clustree:
    input:
        "analyses/{sample}/processed/{sample}.rds"
    output:
        "results/{sample}/clusteringTree/clusteringTree-{sample}.pdf"
    shell:
        "workflow/scripts/scrna-clusteringtree.R --rds {input} --sampleid {wildcards.sample}"
    
rule clustermarkers:
    input:
        "analyses/{sample}/processed/{sample}.rds"
    output:
        "results/{sample}/resolution-{res}/{sample}.positive-markers-forAllClusters.xlsx",
        "results/{sample}/resolution-{res}/{sample}.all-markers-forAllClusters.xlsx"
    shell:
        """
        workflow/scripts/scrna-find-markers.R --rds {input} --resolution {wildcards.res} --sampleid {wildcards.sample}
        """