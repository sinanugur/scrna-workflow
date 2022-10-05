# scRNA sequencing analysis workflow
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) 

Introduction
------------

This is a scRNA-seq analysis pipeline. The pipeline is built in Snakemake and can be run on different platforms and high performance computing (HPC) systems.

Installation
------------

```
git clone https://github.com/sinanugur/scrna-workflow.git
cd scrna-workflow
conda env create -n scrna-workflow -f environment.yml
conda activate scrna-workflow
```

To install R packages:
```
./install_r_packages.sh
```


Quick Start Example
-------------------

The workflow expects Visum 10x samples under data folder in this format:

__"data/{sample}/outs/filtered_feature_bc_matrix.h5"__

This will register the directory name as sample name for later processing.

You can start the pipeline by calling,

```
snakemake -j 20 --configfile config.yaml

```

Do a dry run:
```
snakemake -j 20 -n
```

Show command line arguments:
```
snakemake -j 20 -n -p
```



which will create a 20 threads job.
