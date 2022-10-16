rule create_param_table:
    output:
        tsv="analyses/resolution_detection/{sample}/params.tsv"
    run:
        par_df[["MT","resolution"]].to_csv(output,index=False,sep="\t")