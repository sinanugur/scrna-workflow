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
```

For Apple Silicon (i.e. M1, M2 etc.) architecture, you have to put CONDA_SUBDIR=osx-64 before creating the environment.
```
CONDA_SUBDIR=osx-64 mamba env create --name scrna-workflow --file environment.yml
conda activate scrna-workflow
```


After the environment created and activated successfully, to install all the required R packages, you should run the installation script, this may take some time:
```
bash install_r_packages.sh
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

Now it is time to work on the integrated sample. We can run standard workflow on the integrated object which is always generates at the same location.
```shell
snakemake -j 10 --config  datafolder=analyses_integrated/seurat/integrated.rds option=standard is_integrated_sample=True --rerun-incomplete
```

You may change some of the options or you may provide a config file as well, for example.
```shell
snakemake -j 10 --config  datafolder=analyses_integrated/seurat/integrated.rds option=standard is_integrated_sample=True --configfile=config.yaml --rerun-incomplete
```




