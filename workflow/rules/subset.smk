rule subset_final_rds:
    input:
        rds=files
    output:
        rds="dataoutput/{outname}.rds"
    shell:
        "{cellsnake_path}workflow/scripts/scrna-subset_final_rds.R --rds {input.rds} --output.rds {output.rds} --keywords {keywords} --column {column} --metadata {metadata} {exact}"