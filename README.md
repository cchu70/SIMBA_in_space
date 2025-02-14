# SIMBA_in_space

Create conda environment. Documentation here: https://simba-bio.readthedocs.io/en/latest/Installation.html
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

conda create -n env_simba simba

# in cluster environment
conda create --prefix /data/pinello/SHARED_SOFTWARE/anaconda_latest/envs/cc_envs/env_simba 
conda install python=3.8

# set path for conda to find envs
conda config --add envs_dirs /data/pinello/SHARED_SOFTWARE/anaconda_latest/envs/cc_envs

# install
conda activate env_simba
conda install -c bioconda simba
```

Load environment and install ipykernel for jupyter
```
conda activate env_simba
pip install ipykernel
python -m ipykernel install --user --name env_simba
```

Select `env_simba` as the kernel for running jupyter notebooks.

# Installing R in jupyter lab

- https://docs.ncsa.illinois.edu/systems/delta/en/latest/user_guide/ood/custom-r.html
- https://irkernel.github.io/installation/

1. Open terminal
2. Run `R`
3. `install.packages('IRkernel')`
4. `IRkernel::installspec()`
5. Quit R: `q()`
6. `conda config --add channels r`
7. Open Launcher. `R` should be available

# Installing SpatialPCA on Mac

- Updated to `R version 4.4.2 (2024-10-31)`
- Install `SPARK` from source: `devtools::install_github('xzhoulab/SPARK')`

Issues with installing SPARK on mac. Reported issues for M1 chips: https://github.com/xzhoulab/SPARK/issues/12#issuecomment-2107495029
```
brew install gcc
```

Make files:
```
mkdir -p ~/.R  
nano ~/.R/Makevars
```

From Claude.ai
In `~/.R/Makevars`
```
FC = /usr/local/bin/gfortran
F77 = /usr/local/bin/gfortran
FLIBS = -L/usr/local/Cellar/gcc/14.2.0_1/lib/gcc/current -lgfortran -lquadmath -lm
```

Added following to .bashrc
```
export PATH="/usr/local/bin:$PATH"
```

Add symlinks
```
sudo ln -s /usr/local/Cellar/gcc/14.2.0_1/lib/gcc/current/libgfortran.dylib /usr/local/lib/libgfortran.dylib
sudo ln -s /usr/local/Cellar/gcc/14.2.0_1/lib/gcc/current/libquadmath.dylib /usr/local/lib/libquadmath.dylib
```

In an R session:
```
R
> install.packages("devtools")
> library(devtools)
> devtools::install_github("xzhoulab/SPARK") # install separately
> install_github("shangll123/SpatialPCA")
> remotes::install_github("mojaveazure/seurat-disk") # to read h5ad files
```

## Install on cluster

Create R conda environment
```
conda create --prefix /data/pinello/SHARED_SOFTWARE/anaconda_latest/envs/cc_envs/R_env
conda activate R_env
conda install -c conda-forge r-base
conda install -c conda-forge r-essentials
conda install conda-forge::r-rspectra

conda install python==3.9 # requirement for r-seurat?
conda install bioconda::r-seurat

#conda install conda-forge::r-umap
#conda install bioconda::bioconductor-splatter
```

Running into issues installing r-seurat. Claude.ai 3.5 Sonnet:
```
conda install -c conda-forge r-matrix
conda install -c conda-forge r-rcpp
conda install -c bioconda r-uwot
conda install -c bioconda r-seurat
```

Tried remaking from scratch:
```
conda create -n cc_simba r-base=4.3 python=3.10 -c conda-forge
# Error/Missing dependencies


R
> install.packages(c("SeuratObject", "cowplot", "fitdistrplus", "ggplot2", "ggrepel", "ggridges", "igraph", "irlba", "leidenbase", "MASS", "Matrix", "miniUI", "patchwork", "plotly", "png", "reticulate", "RSpectra", "scattermore", "sctransform", "shiny", "spatstat.explore", "spatstat.geom", "uwot", "RcppEigen"))

Error: 
ERROR: dependencies ‘spatstat.data’, ‘spatstat.geom’, ‘spatstat.random’, ‘spatstat.sparse’, ‘Matrix’ are not available for package ‘spatstat.explore’

conda install conda-forge::r-matrix
conda install conda-forge::r-spatstat.geom
conda install conda-forge::r-spatstat.random conda-forge::r-spatstat.sparse conda-forge::r-spatstat.data
conda install conda-forge::r-spatstat.explore
conda install conda-forge::r-ggplot2
conda install conda-forge::r-rcppeigen conda-forge::r-uwot
conda install conda-forge::r-shiny conda-forge::r-sctransform conda-forge::r-scattermore conda-forge::r-rspectra conda-forge::r-reticulate
conda install conda-forge::r-png conda-forge::r-plotly conda-forge::r-patchwork conda-forge::r-miniui conda-forge::r-mass bioconda::r-leidenbase
conda install r::r-irlba
conda install conda-forge::r-igraph
conda install conda-forge::r-ggridges conda-forge::r-ggrepel conda-forge::r-ggplot2 conda-forge::r-fitdistrplus

conda install conda-forge::r-cowplot conda-forge::r-seuratobject
conda install -c bioconda r-seurat

conda install conda-forge::r-devtools

# Error loading devtools library in R shell
conda install conda-forge::r-fastmap=1.2.0
conda install conda-forge::r-promises=1.3.2
conda install bioconda::bioconductor-splatter

# more requirements for SpatialPCA tutorial
conda install bioconda::bioconductor-bluster

# simba
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install bioconda::simba --no-deps # taking a long time. too many dependencies?
# try to make a separate environment for this
conda deactivate # in lab_py38 env
conda create -n cc_simba_only simba # took a long time

# cloning existing environment
conda create --name cc_simba_only --clone jy_simba
conda activate cc_simba_only
(cc_simba_only): python -m ipykernel install --user --name cc_simba_only
(cc_simba_only): pip install --upgrade matplotlib # upgrade to 3.10 for scanpy plotting compatibility

# running scanpy leiden clustering
(cc_simba_only): pip3 install igraph
(cc_simba_only): pip3 install leidenalg squidpy
```

Install SpatialPCA
```
R
> install.packages("devtools")
> library(devtools)
> install_github("shangll123/SpatialPCA")

# for scripting
> install.packages("argparse") 
> install.packages("optparse")
```

Jupyter
```
conda install ipykernel
python -m ipykernel install --user --name cc_simba --display-name cc_simba
pip install scanpy seaborn matplotlib pandas numpy
```

Follow tutorial: https://lulushang.org/SpatialPCA_Tutorial/DLPFC.html

# Setting up conda environment on cluster

Make the `yml` file:
```
conda env export --no-builds > environment.yml
```

Removed `prefix:` line in the `environment.yml` file, and line with `- libgfortran=5.0.0`
```
conda env create --prefix /data/pinello/SHARED_SOFTWARE/anaconda_latest/envs/cc_envs/env_simba --file=environment.yml
```

Or, just remake:
```
conda env create --prefix /data/pinello/SHARED_SOFTWARE/anaconda_latest/envs/cc_envs/env_simba simba
```

## Using mamba to make simba environment

In ml007
```
../../SHARED_SOFTWARE/anaconda_latest/bin/mamba create -n cc_mamba_simba bioconda::simba

conda activate cc_mamba_simba
../../SHARED_SOFTWARE/anaconda_latest/bin/mamba install squidpy
```

# Running benchmarking analysis

Tasks:
- Embedding tools:
    - PCA
    - SpatialPCA
    - SIMBA RNAseq only
    - SIMBA spatial info only
- Performance metrics
    - Silouhette score
    - ARI/NMI/LISI
 
Data: human DLPFC