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
        input_function
    output:
        rds="analyses/raw/{sample}/" + "MT~{MT}/resolution~{resolution}.rds",
        before="results/{sample}/" + "MT~{MT}/resolution~{resolution}" + "/technicals/before-qc-trimming-violinplot.pdf",
        after="results/{sample}/" + "MT~{MT}/resolution~{resolution}" + "/technicals/after-qc-trimming-violinplot.pdf"
    shell:
        """workflow/scripts/scrna-read-qc.R --data.dir {input} --output.rds {output.rds} --sampleid {wildcards.sample} --percent.mt {wildcards.MT} --min.features {min_features} --min.cells {min_cells} --before.violin.plot {output.before} --after.violin.plot {output.after}"""

