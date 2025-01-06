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
