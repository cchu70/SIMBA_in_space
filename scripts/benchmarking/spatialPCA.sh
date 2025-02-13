#!/bin/bash

# Run from scripts/. directory


conda activate cc_simba # R environment
Rscript benchmarking/spatialPCA.r --adata_dir ../data/human_DLPFC_RData_Shang/ --workdir ../results/00/spatialpca  #--copy_rdata_only
conda deactivate


