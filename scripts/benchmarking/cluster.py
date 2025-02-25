import os
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
import tqdm

# approximate original figure in http://spatial.libd.org/spatialLIBD/
palette_celltype={'L1':'#eb34a8',
                  'L2':'#3486eb',
                  'L3':'#34eb5b',
                  'L4':"#ae34eb",
                  'L5':'#ebdb34',
                  'L6':'#eb9234',
                  'WM':'#000000'}

palette_entity_anno = palette_celltype.copy()
palette_entity_anno['gene'] = "lightgray"

import simba as si

si.settings.set_figure_params(
    dpi=80,
    style='white',
    fig_size=[5,5],
    rc={'image.cmap': 'viridis'}
)

from matplotlib_inline.backend_inline import set_matplotlib_formats
set_matplotlib_formats('retina')
import matplotlib.pyplot as plt

from sklearn.metrics import silhouette_score, adjusted_rand_score, normalized_mutual_info_score

def get_pairwise_distance(adata_matrix):
    # adata_matrix: adata.X, or adata.obsm['X_PCA'], etc.
    diff = adata_matrix[:, np.newaxis, :] - adata_matrix[np.newaxis, :, :]
    distances = np.sqrt(np.sum(diff**2, axis=2) )
    return distances 

def run_leiden(
    adata, # adata_CG
    true_label_col='spatialLIBD'
):
    # some cells have nans in the spatialLIBD column in spatial PCA
    nan_true_labels = adata.obs[true_label_col].isna()
    idx = adata.obs[~nan_true_labels].index.tolist()
    adata = adata[idx].copy()

    # Silhouette score
    dists = get_pairwise_distance(adata.X) # on the data in X
    sil = silhouette_score(X=dists, labels=adata.obs[true_label_col])

    # perform leiden clustering
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata)

    # scores
    ari = adjusted_rand_score(adata.obs[true_label_col], adata.obs['leiden'])
    nmi = normalized_mutual_info_score(adata.obs[true_label_col], adata.obs['leiden'])

    return adata, sil, ari, nmi


def pca_leiden(
    adata, # adata_CG
    true_label_col='spatialLIBD'
):
    # some cells have nans in the spatialLIBD column in spatial PCA
    nan_true_labels = adata.obs[true_label_col].isna()
    idx = adata.obs[~nan_true_labels].index.tolist()
    adata = adata[idx].copy()

    # Silhouette score
    pca_dists = get_pairwise_distance(adata.obsm['X_pca']) 
    pca_sil = silhouette_score(X=pca_dists, labels=adata.obs[true_label_col])

    # Copy PCA to X
    pc_adata = ad.AnnData(
        X=adata.obsm['X_pca'],
        obs=adata.obs,
    )
    pc_adata.obsm['X_pca'] = adata.obsm['X_pca'].copy()

    # perform leiden clustering
    sc.pp.neighbors(pc_adata)
    sc.tl.leiden(pc_adata)

    # scores
    ari = adjusted_rand_score(pc_adata.obs[true_label_col], pc_adata.obs['leiden'])
    nmi = normalized_mutual_info_score(pc_adata.obs[true_label_col], pc_adata.obs['leiden'])

    return pc_adata, pca_sil, ari, nmi

def run_leiden_clustering(
    adata_output_df,
    performance_output_fn,
    path_col='run_pca',
    cell_embedding_adata_fn="adata_CG.h5ad",
    true_label_col='spatialLIBD',
    version='normal', # 'PCA'
):
    performance_df = pd.DataFrame(index=adata_output_df.index, columns=['silhoutte', 'ARI', 'NMI', 'leiden_adata_fn'])
    for sample, adata_dir in tqdm.tqdm(adata_output_df[path_col].items(), total=adata_output_df.shape[0]):
        adata_fn = f"{adata_dir}/{cell_embedding_adata_fn}"
        leiden_adata_fn = f"{adata_dir}/leiden.{cell_embedding_adata_fn}"

        try:
            adata = sc.read_h5ad(adata_fn)
        except:
            print(f"{adata_fn} does not exist. Skipping")
            continue

        if version == 'PCA':
            leiden_adata, pca_sil, ari, nmi = pca_leiden(adata, true_label_col=true_label_col)
        else:
            leiden_adata, pca_sil, ari, nmi = run_leiden(adata, true_label_col=true_label_col)
        
        leiden_adata.write(leiden_adata_fn)

        performance_df.loc[sample, 'silhoutte'] = pca_sil
        performance_df.loc[sample, 'ARI'] = ari
        performance_df.loc[sample, 'NMI'] = nmi
        performance_df.loc[sample, 'leiden_adata_fn'] = leiden_adata_fn

    performance_df.to_csv(performance_output_fn, sep='\t')


def run_walktrap(
    adata, # adata_CG
    true_label_col='spatialLIBD'
):
     # some cells have nans in the spatialLIBD column in spatial PCA
    nan_true_labels = adata.obs[true_label_col].isna()
    idx = adata.obs[~nan_true_labels].index.tolist()
    adata = adata[idx].copy()

    # perform leiden clustering
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata)

    # scores
    ari = adjusted_rand_score(adata.obs[true_label_col], adata.obs['leiden'])
    nmi = normalized_mutual_info_score(adata.obs[true_label_col], adata.obs['leiden'])

    return adata, ari, nmi
