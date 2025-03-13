from benchmarking.niche import run_niche_comparison

sb_global_performance_fn, sb_global_figures_fn = run_niche_comparison(
    results_tsv_path='../results/00/simba_spatial_and_BANKSY/mask_e500_l0.8_b10/run_simba_spatial_and_rna_BANKSY.output.tsv', 
    cell_embedding_adata_fn="adata_C.h5ad",
    adata_spatial_dir="/data/pinello/PROJECTS/2025-01-31_CC_Spatial_SIMBA/SIMBA_in_space/data/human_DLPFC/",
    col_name='run_simba_spatial_and_rna_BANKSY', # column name with path to sample results
)


# from scipy.spatial import cKDTree
# from scipy.stats import spearmanr, pearsonr
# from sklearn.metrics.pairwise import cosine_similarity
# from sklearn.manifold import trustworthiness
# import seaborn as sns
# import tqdm
# import anndata as ad
# import numpy as np
# import pandas as pd
# from scipy.sparse import csr_matrix
# import matplotlib.pyplot as plt
# import os
# import scanpy as sc

