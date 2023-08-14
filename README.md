# scRNA sequencing analysis workflow
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) 

Introduction
------------

Cellsnake can be run directly using the snakemake workflow. We recommend the wrapper but the snakemake workflow give more control in some use cases.

The main cellsnake repo is here : https://github.com/sinanugur/cellsnake


Installation
------------

You may pull the workflow from the GitHub repo and create a clean environment. Mamba installation is highly recommended.

```
conda install mamba -c conda-forge # to install Mamba

git clone https://github.com/sinanugur/scrna-workflow.git
cd scrna-workflow
mamba env create --name scrna-workflow --file environment.yml
conda activate scrna-workflow

mamba env create --name cellsnake_testing --file environment.yml
```

For Apple Silicon (i.e. M1, M2 etc.) architecture, you have to put CONDA_SUBDIR=osx-64 before creating the environment.
```
CONDA_SUBDIR=osx-64 mamba env create --name scrna-workflow --file environment.yml
```


After the environent created and activated, to install required R packages:
```
./install_r_packages.sh
```


Quick Start Example
-------------------

You can start a minimal run by calling, sample runs are expected in data folder.

```shell
snakemake -j 10 --config datafolder=data option=minimal
```

Then we can run integration.
```shell
snakemake -j 10 --config option=integration
```

Now it is time to work on the integrated sample. We can run full advanced run on the integrated object which is always generates at the same location.
```shell
snakemake -j 10 --config  datafolder=analyses_integrated/seurat/integrated.rds resolution=0.3 option=advanced is_integrated_sample=True --rerun-incomplete
```





Do a dry run:
```
snakemake -j 5 -n
```

Show command line arguments and dry run:
```
snakemake -j 5 -n -p
```
