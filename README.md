# SIMBA_in_space

Create conda environment. Documentation here: https://simba-bio.readthedocs.io/en/latest/Installation.html
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

conda create -n env_simba simba
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

# Installing SpatialPCA

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

In `~/.R/Makevars`
```
FC = /usr/local/bin/gfortran
F77 = /usr/local/bin/gfortran
#FLIBS = -L/usr/local/Cellar/gcc/14.2.0_1 -lgfortran -lquadmath -lm
FLIBS = /usr/local/Cellar/gcc/14.2.0_1 -lgfortran -lquadmath -lm
```

Added following to .bashrc
```
export PATH="/usr/local/bin:$PATH"
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