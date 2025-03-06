#!/bin/bash

# Run from scripts/. directory
# activate simba conda environment (cc_mamba_simba)

# <method>_r<radius/gamma>_s<pearson scaling>_w<CC edge type training weight>

# python benchmarking/simba_spatial_and_rna.py \
#     --workdir ../results/00/simba_spatial_and_rna \
#     --adata_dir ../data/human_DLPFC \
#     --spatial_method rbf \
#     --e 5 \
#     --CC_training_weight 1.0 \
#     --pearson_corr_scale 

# Rerun one sample

export PATH="/data/pinello/SHARED_SOFTWARE/anaconda_latest/envs/cc_envs/cc_mamba_simba/bin:$PATH"

python benchmarking/simba_spatial_and_rna.py \
    --workdir ../results/00/simba_spatial_and_rna \
    --adata_fns ../data/human_DLPFC/151507.h5ad \
    --spatial_method rbf \
    --e 5 \
    --CC_training_weight 1.0 \
    --pearson_corr_scale 