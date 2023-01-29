
def kraken2_input_function(wildcards):
    if os.path.isfile(datafolder + "/" + wildcards.sample + "/possorted_genome_bam.bam"):
        return(datafolder +  "/" + wildcards.sample + "/possorted_genome_bam.bam")
    elif os.path.isfile(datafolder + "/" + wildcards.sample + "/outs/possorted_genome_bam.bam"):
        print("burdayim")
        return(datafolder + "/" + wildcards.sample + "/outs/possorted_genome_bam.bam")
    else:
        return


rule run_kraken:
    input:
        bam=kraken2_input_function
    output:
        analyses_folder + "/kraken/" + f"{paramspace.wildcard_pattern}" + "/{sample}/counts/matrix.mtx",
        analyses_folder + "/kraken/" + f"{paramspace.wildcard_pattern}" + "/{sample}/counts/hierarchy.txt",
        directory(analyses_folder + "/kraken/" + f"{paramspace.wildcard_pattern}" + "/{sample}/")
    threads: 5
 
    shell:
        "workflow/mg2sc/src/scMeG-kraken.py --input {input.bam} --outdir {ouput[2]} --DBpath {kraken_db_folder} --threads {threads} --prefix {wildcards.sample}"


rule collapse_kraken:
    input:
        outdir=analyses_folder + "/kraken/" + f"{paramspace.wildcard_pattern}" + "/{sample}/",
        hierarchy=analyses_folder + "/kraken/" + f"{paramspace.wildcard_pattern}" + "/{sample}/counts/hierarchy.txt"
    output:
        analyses_folder + "/kraken/" + f"{paramspace.wildcard_pattern}" + "/{sample}/{taxa}.h5ad"

    shell:
        "workflow/scripts/scrna-kraken2-collapse.py {input.outdir} {input.hierarchy} {wilcards.taxa} {output}"

rule convert_to_seurat:
    input:
        analyses_folder + "/kraken/" + f"{paramspace.wildcard_pattern}" + "/{sample}/{taxa}.h5ad"
    output:
        analyses_folder + "/kraken/" + f"{paramspace.wildcard_pattern}" + "/{sample}/{taxa}.h5seurat"
    shell:
        """Rscript -e 'SeuratDisk::Convert({input},dest = "h5seurat", overwrite = TRUE)'"""