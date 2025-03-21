#!/bin/bash

# Run from scripts/. directory
# activate simba conda environment (cc_mamba_simba)

python benchmarking/simba_spatial_only.py \
    --workdir ../results/00/simba_spatial_only \
    --adata_dir ../data/human_DLPFC \
    --spatial_method rbf \
    --e 5