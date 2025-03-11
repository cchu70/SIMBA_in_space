#!/bin/bash

# Run from scripts/. directory
# activate simba conda environment (cc_mamba_simba)

# <method>_g<grid size scaling>_w<CR edge type training weight>

python benchmarking/simba_spatial_and_rna_rankK.py \
    --workdir ../results/00/simba_spatial_and_rna_rankK \
    --adata_dir ../data/human_DLPFC \
    --lambda_val 0.8 \
    --k 6 \
    --n_bins 20 
