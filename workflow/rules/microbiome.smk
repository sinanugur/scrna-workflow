
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
        matrix=analyses_folder + "/kraken/" + f"{paramspace.wildcard_pattern}" + "/{sample}/counts/matrix.mtx",
        hierarchy=analyses_folder + "/kraken/" + f"{paramspace.wildcard_pattern}" + "/{sample}/counts/hierarchy.txt",
        outdir=directory(analyses_folder + "/kraken/" + f"{paramspace.wildcard_pattern}" + "/{sample}/"),
        unmapped=temp(analyses_folder + "/kraken/" + f"{paramspace.wildcard_pattern}" + "/{sample}/{sample}" + "_unmapped.bam"),
        fq=temp(analyses_folder + "/kraken/" + f"{paramspace.wildcard_pattern}" + "/{sample}/{sample}" + "_unmapped.fq"),
        kr=temp(analyses_folder + "/kraken/" + f"{paramspace.wildcard_pattern}" + "/{sample}/{sample}" + "_output.kraken")
    threads: 10
    shell:
        """
        rm -r {output.outdir}/counts
        workflow/mg2sc/src/scMeG-kraken.py --input {input.bam} --outdir {output.outdir} --DBpath {kraken_db_folder} --threads {threads} --prefix {wildcards.sample}
        """


rule collapse_kraken:
    input:
        outdir=analyses_folder + "/kraken/" + f"{paramspace.wildcard_pattern}" + "/{sample}",
        hierarchy=analyses_folder + "/kraken/" + f"{paramspace.wildcard_pattern}" + "/{sample}/counts/hierarchy.txt"
    output:
        analyses_folder + "/kraken/" + f"{paramspace.wildcard_pattern}" + "/{sample}/{taxa}.h5ad"

    shell:
        "workflow/scripts/scrna-kraken2-collapse.py {input.outdir}/counts {input.hierarchy} {wildcards.taxa} {output}"

rule convert_to_seurat:
    input:
        analyses_folder + "/kraken/" + f"{paramspace.wildcard_pattern}" + "/{sample}/{taxa}.h5ad"
    output:
        analyses_folder + "/kraken/" + f"{paramspace.wildcard_pattern}" + "/{sample}/{taxa}.h5seurat"
    shell:
        """Rscript -e 'SeuratDisk::Convert("{input}",dest = "h5seurat", overwrite = TRUE)'"""

rule parse_h5seurat:
    input:
        analyses_folder + "/kraken/" + f"{paramspace.wildcard_pattern}" + "/{sample}/{taxa}.h5seurat"
    output:
        microbiome_rds=analyses_folder + "/kraken/" + f"{paramspace.wildcard_pattern}" + "/{sample}/microbiome-full-{taxa}-level.rds",
        plot=results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/microbiome/microbiome-full-{taxa}-level.pdf"
    shell:
        "workflow/scripts/scrna-kraken2-data-parser.R --h5seurat {input} --output.rds {output.microbiome_rds} --output.plot {output.plot} --sampleid {wildcards.sample} --taxa {wildcards.taxa} --min.features {microbiome_min_features} --min.cells {microbiome_min_cells}"


rule dimplot_for_microbiome:
    input:
        rds=analyses_folder + "/processed/" + f"{paramspace.wildcard_pattern}" + "/{sample}.rds",
        microbiome_rds=analyses_folder + "/kraken/" + f"{paramspace.wildcard_pattern}" + "/{sample}/microbiome-full-{taxa}-level.rds"
    output:
        results_folder + "/{sample}/" + f"{paramspace.wildcard_pattern}" + "/microbiome/dimplot-{taxa}-{reduction}.pdf"
    shell:
        "workflow/scripts/scrna-microbiome-dimplot.R --rds {input.rds} --microbiome.rds {input.microbiome_rds} --output.plot {output} --reduction.type {wildcards.reduction} --taxa {wildcards.taxa}"


