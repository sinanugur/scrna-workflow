rule integration_with_seurat:
    input:
        integration_files
    output:
        "analyses_integrated/seurat/" + integration_id + ".rds"
    shell:
        """
        {cellsnake_path}workflow/scripts/scrna-seurat-integration.R --rds "{input}" --sampleid {integration_id} --output.rds {output} --reduction {reduction} --dims {dims}
        """

