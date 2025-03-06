#!/bin/bash

# Run from scripts/. directory
# activate simba conda environment (cc_mamba_simba)

# <method>_g<grid size scaling>_w<CR edge type training weight>

export PATH="/data/pinello/SHARED_SOFTWARE/anaconda_latest/envs/cc_envs/cc_mamba_simba/bin:$PATH"

python benchmarking/simba_spatial_and_rna_BANKSY.py \
    --workdir ../results/00/simba_spatial_and_BANKSY \
    --adata_fns ../data/human_DLPFC/151507.h5ad \
    --lambda_val 0.8 \
    --spatial_method mask \
    --e 500 \
    --n_bins 10 \