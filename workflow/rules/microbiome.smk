
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
        outdir=directory(analyses_folder + "/kraken/" + f"{paramspace.wildcard_pattern}" + "/{sample}/")
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