#!/bin/bash

# Run from scripts/. directory
# activate simba conda environment

python benchmarking/simba_spatial_only.py \
    --workdir ../results/00/simba_spatial_only \
    --adata_dir ../data/human_DLPFC \
    --e 3 \
    --rerun