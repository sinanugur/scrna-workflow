rule subset_final_rds:
    input:
        rds=files
    output:
        rds="dataoutput/{subset}.rds"
    shell:
        "{cellsnake_path}workflow/scripts/scrna-subset_final_rds.R --rds {input.rds} --output.rds {output.rds} --keywords {keywords} --column {subset_column} --metadata {metadata} {exact}"