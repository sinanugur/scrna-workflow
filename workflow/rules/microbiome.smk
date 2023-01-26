
def kraken2_input_function(wildcards):
    if os.path.isfile(datafolder):
        return(datafolder)
    else:
        if os.path.isfile(datafolder + "/" + wildcards.sample + "/possorted_genome_bam.bam"):
            return(datafolder +  "/" + wildcards.sample + "/")
        elif os.path.isfile(datafolder + "/" + wildcards.sample + "/outs/possorted_genome_bam.bam"):
            return(datafolder + "/" + wildcards.sample + "/outs")




