rule auto_detect_resolution:
    input:
        input_function
    output:
        tsv="analyses/resolution_detection/{sample}/" + f"{paramspace_mt.wildcard_pattern}.tsv"
    params:
        paramaters=paramspace_mt.instance,
    threads:5

    shell:
        "{cellsnake_path}workflow/scripts/scrna-cluster-resolution-detection.R --data.dir {input} --sampleid {wildcards.sample} --percent.mt {params.paramaters[MT]} --min.features {min_features} --min.cells {min_cells} --output.tsv {output.tsv} --cpu {threads}"


rule create_param_table:
    input:
        lambda wildcards: expand("analyses/resolution_detection/" + wildcards.sample  + "/{params}.tsv",params=list(paramspace_mt.instance_patterns))
    output:
        tsv="analyses/resolution_detection/{sample}/params.tsv"
    run:
        with open(output.tsv,"w") as o:
            par_df =  pd.DataFrame()
            for i in input:
                par_df=pd.concat([par_df,pd.read_table(i)])
            
            par_df[["MT","resolution"]].to_csv(o,index=False,sep="\t")


            