#!/bin/bash

# Run from scripts/. directory
# activate simba conda environment (cc_mamba_simba)

# <method>_r<radius/gamma>_s<pearson scaling>_w<CC edge type training weight>

python benchmarking/simba_spatial_and_rna.py \
    --workdir ../results/00/simba_spatial_and_rna \
    --adata_dir ../data/human_DLPFC \
    --spatial_method mask \
    --e 500 \
    --CC_training_weight 1.0 \
    --pearson_corr_scale \
    --rerun