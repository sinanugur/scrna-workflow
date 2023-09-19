
def kraken2_input_function(wildcards):
    if os.path.isfile(datafolder + "/" + wildcards.sample + "/possorted_genome_bam.bam"):
        return(datafolder +  "/" + wildcards.sample + "/possorted_genome_bam.bam")
    elif os.path.isfile(datafolder + "/" + wildcards.sample + "/outs/possorted_genome_bam.bam"):
        return(datafolder + "/" + wildcards.sample + "/outs/possorted_genome_bam.bam")
    else:
        return(datafolder + "/" + wildcards.sample + "/outs/possorted_genome_bam.bam")





rule run_kraken:
    input:
        bam=kraken2_input_function
    output:
        matrix=analyses_folder + "/kraken/" + "{confidence}_{min_hit_groups}/" + f"{paramspace.wildcard_pattern}" + "/{sample}/counts/matrix.mtx",
        hierarchy=analyses_folder + "/kraken/" + "{confidence}_{min_hit_groups}/" + f"{paramspace.wildcard_pattern}" + "/{sample}/counts/hierarchy.txt",
        outdir=directory(analyses_folder + "/kraken/" + "{confidence}_{min_hit_groups}/"  + f"{paramspace.wildcard_pattern}" + "/{sample}/"),
        unmapped=temp(analyses_folder + "/kraken/" + "{confidence}_{min_hit_groups}/" + f"{paramspace.wildcard_pattern}" + "/{sample}/{sample}" + "_unmapped.bam") if kraken_extra_files is False else analyses_folder + "/kraken/" + "{confidence}_{min_hit_groups}/" + f"{paramspace.wildcard_pattern}" + "/{sample}/{sample}" + "_unmapped.bam",
        fq=temp(analyses_folder + "/kraken/" + "{confidence}_{min_hit_groups}/" + f"{paramspace.wildcard_pattern}" + "/{sample}/{sample}" + "_unmapped.fq") if kraken_extra_files is False else analyses_folder + "/kraken/" + "{confidence}_{min_hit_groups}/" + f"{paramspace.wildcard_pattern}" + "/{sample}/{sample}" + "_unmapped.fq",
        kr=analyses_folder + "/kraken/" + "{confidence}_{min_hit_groups}/"  + f"{paramspace.wildcard_pattern}" + "/{sample}/{sample}" + "_output.kraken" if kraken_extra_files is True else temp(analyses_folder + "/kraken/" + "{confidence}_{min_hit_groups}/"  + f"{paramspace.wildcard_pattern}" + "/{sample}/{sample}" + "_output.kraken"),
        classified=analyses_folder + "/kraken/" + "{confidence}_{min_hit_groups}/"  + f"{paramspace.wildcard_pattern}" + "/{sample}/{sample}" + "_classified_sequences.txt" if kraken_extra_files is True else temp(analyses_folder + "/kraken/" + "{confidence}_{min_hit_groups}/"  + f"{paramspace.wildcard_pattern}" + "/{sample}/{sample}" + "_classified_sequences.txt")
    threads: 10
    params:
        bowtie= "--bowtie " + bowtie_database_prefix  if bowtie_database_prefix is not None else ""
        predb= "--predb " + prekraken_db_folder  if prekraken_db_folder is not None else ""
    run:
        shell("""
            rm -r {output.outdir}/counts;
            {cellsnake_path}workflow/mg2sc/src/scMeG-kraken.py --input {input.bam} {params.bowtie} {params.predb} --outdir {output.outdir} --DBpath {kraken_db_folder} --threads {threads} --minimum-hit-groups {min_hit_groups} --confidence {confidence} --complexity {complexity} --prefix {wildcards.sample}
            rm -f {output.unmapped}.tmp
            rm -f {output.fq}.tmp
            rm -f {output.fq}.id
            rm -f {output.fq}.exclude
            rm -f {output.kr}.tmp
            """)
        #if kraken_extra_files is False:
        #    shell("rm -r {output.kr} {output.classified}")


rule collapse_kraken:
    input:
        outdir=analyses_folder + "/kraken/" + "{confidence}_{min_hit_groups}/"  + f"{paramspace.wildcard_pattern}" + "/{sample}",
        hierarchy=analyses_folder + "/kraken/" + "{confidence}_{min_hit_groups}/"  + f"{paramspace.wildcard_pattern}" + "/{sample}/counts/hierarchy.txt"
    output:
        analyses_folder + "/kraken/" + "{confidence}_{min_hit_groups}/"  + f"{paramspace.wildcard_pattern}" + "/{sample}/{taxa}.h5ad"

    shell:
        "{cellsnake_path}workflow/scripts/scrna-kraken2-collapse.py {input.outdir}/counts {input.hierarchy} {wildcards.taxa} {output}"

rule convert_to_seurat:
    input:
        analyses_folder + "/kraken/" + "{confidence}_{min_hit_groups}/"  + f"{paramspace.wildcard_pattern}" + "/{sample}/{taxa}.h5ad"
    output:
        analyses_folder + "/kraken/" + "{confidence}_{min_hit_groups}/"  + f"{paramspace.wildcard_pattern}" + "/{sample}/{taxa}.h5seurat"
    shell:
        """Rscript -e 'SeuratDisk::Convert("{input}",dest = "h5seurat", overwrite = TRUE)'"""

rule parse_h5seurat:
    input:
        analyses_folder + "/kraken/" + "{confidence}_{min_hit_groups}/"  + f"{paramspace.wildcard_pattern}" + "/{sample}/{taxa}.h5seurat"
    output:
        microbiome_rds=analyses_folder + "/kraken/" + "{confidence}_{min_hit_groups}/" + f"{paramspace.wildcard_pattern}" + "/{sample}/microbiome-full-{taxa}-level.rds",
        plot=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/microbiome/" + "{confidence}_{min_hit_groups}/" + "plot_microbiome_barplot-{taxa}-level.pdf"
    shell:
        "{cellsnake_path}workflow/scripts/scrna-kraken2-data-parser.R --h5seurat {input} --output.rds {output.microbiome_rds} --output.plot {output.plot} --sampleid {wildcards.sample} --taxa {wildcards.taxa} --min.features {microbiome_min_features} --min.cells {microbiome_min_cells}"


rule dimplot_for_microbiome:
    input:
        rds=analyses_folder + "/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds",
        microbiome_rds=analyses_folder + "/kraken/" + "{confidence}_{min_hit_groups}/"  + f"{paramspace.wildcard_pattern}" + "/{sample}/microbiome-full-{taxa}-level.rds" if is_integrated_sample is False else "analyses_integrated/seurat/" + "{confidence}_{min_hit_groups}/"  + integration_id + "-{taxa}.rds"
    output:
        dimplot=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/microbiome/" + "{confidence}_{min_hit_groups}/" + "plot_microbiome_dimplot-{taxa}-{reduction}.pdf",
        tplot=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/microbiome/" + "{confidence}_{min_hit_groups}/" + "plot_microbiome_total_dimplot-{taxa}-{reduction}.pdf"
    shell:
        "{cellsnake_path}workflow/scripts/scrna-microbiome-dimplot.R --rds {input.rds} --microbiome.rds {input.microbiome_rds} --dimplot {output.dimplot} --tplot {output.tplot} --reduction.type {wildcards.reduction} --taxa {wildcards.taxa}"


rule table_for_microbiome:
    input:
        rds=analyses_folder + "/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds",
        microbiome_rds=analyses_folder + "/kraken/" + "{confidence}_{min_hit_groups}/"  + f"{paramspace.wildcard_pattern}" + "/{sample}/microbiome-full-{taxa}-level.rds" if is_integrated_sample is False else "analyses_integrated/seurat/" + "{confidence}_{min_hit_groups}/"  + integration_id + "-{taxa}.rds"

    output:
        #dimplot=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/microbiome/plot_significance-{taxa}-{i}.pdf",
        xlsx=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/microbiome/" + "{confidence}_{min_hit_groups}/" + "table_microbiome-{taxa}-{i}.xlsx"
    shell:
        "{cellsnake_path}workflow/scripts/scrna-microbiome-sigplot.R --rds {input.rds} --microbiome.rds {input.microbiome_rds} --taxa {wildcards.taxa}  --sigtable {output.xlsx} --idents {wildcards.i}"



rule combine_microbiome_files_for_later:
    input:
        lambda wildcards : expand([analyses_folder + "/kraken/" + wildcards.confidence + "_" + wildcards.min_hit_groups + "/"  + "{params}" + "/" + s + "/microbiome-full-" + wildcards.taxa + "-level.rds" for s in files],params=list(paramspace.instance_patterns))
    output:
        "analyses_integrated/seurat/" + "{confidence}_{min_hit_groups}/" + integration_id + "-{taxa}.rds"
    shell:
        """{cellsnake_path}workflow/scripts/scrna-combine-microbiome-rds.R --rds "{input}" --sampleid {integration_id} --output.rds {output}"""
