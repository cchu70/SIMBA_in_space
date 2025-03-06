#!/bin/bash

# Run from scripts/. directory
# activate simba conda environment (cc_mamba_simba)

# <method>_g<grid size scaling>_w<CR edge type training weight>
# --adata_dir ../data/human_DLPFC \

export PATH="/data/pinello/SHARED_SOFTWARE/anaconda_latest/envs/cc_envs/cc_mamba_simba/bin:$PATH"

python benchmarking/simba_spatial_and_rna_regions_entities.py \
    --workdir ../results/00/simba_spatial_and_rna_region_entities \
    --adata_fns ../data/human_DLPFC/151507.h5ad \
    --grid_size 500 \
    --CR_training_weight 1.0 \
    --rerun