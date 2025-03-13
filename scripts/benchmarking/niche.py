from scipy.spatial import cKDTree
from scipy.stats import spearmanr, pearsonr
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.manifold import trustworthiness
import seaborn as sns
import tqdm
import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import os
import scanpy as sc

def get_k_neighbors(adata_CG):
    
    coords = adata_CG.obsm['spatial']
    tree = cKDTree(coords)

    k=6
    # Query the k+1 nearest neighbors (to exclude the point itself)
    distances, indices = tree.query(coords, k=k+1)
    # Exclude self (first column)
    indices = indices[:, 1:]
    indices_df = pd.DataFrame(indices, index=adata_CG.obs.index, columns=np.arange(1, k + 1))
    return indices_df

def get_label_distr(adata_CG, label):
    indices_df = get_k_neighbors(adata_CG)

    cell_type_distr_df = pd.DataFrame(index=adata_CG.obs.index, columns=adata_CG.obs[label].unique())
    for barcode in cell_type_distr_df.index:
        k_neighbors = indices_df.loc[barcode].to_numpy()
        celltype_counts = adata_CG.obs.loc[adata_CG.obs.index[k_neighbors]][label].value_counts()
        cell_type_distr_df.loc[barcode, celltype_counts.index] = celltype_counts
    return cell_type_distr_df

def get_cossim(df):
    full_cossim = cosine_similarity(df)
    cossim_flat = full_cossim[np.triu_indices_from(full_cossim, k=1)]  # remove diagonal
    return cossim_flat

def get_cossim_corr(adata_CG, adata_C, label='spatialLIBD', trust_n_neighbors=6):
    cell_type_distr_df = get_label_distr(adata_CG, label=label)
    cell_type_distr_cossim_flat = get_cossim(cell_type_distr_df)
    embed_cossim_flat = get_cossim(adata_C[cell_type_distr_df.index, :].X)

    spearman_corr, _ = spearmanr(cell_type_distr_cossim_flat, embed_cossim_flat)
    pearson_corr, _ = pearsonr(cell_type_distr_cossim_flat, embed_cossim_flat)
    # Assume features and embeddings are standardized (or not, depending on scale)
    trust = trustworthiness(cell_type_distr_df.to_numpy(), adata_C[cell_type_distr_df.index, :].X, n_neighbors=trust_n_neighbors)

    return dict(cell_type_distr_df=cell_type_distr_df, cell_type_distr_cossim_flat=cell_type_distr_cossim_flat, embed_cossim_flat=embed_cossim_flat), \
        dict(spearman_corr=spearman_corr, pearson_corr=pearson_corr, trust=trust)

def run_niche_comparison(
    results_tsv_path, 
    cell_embedding_adata_fn="adata_C.h5ad",
    adata_spatial_dir="/data/pinello/PROJECTS/2025-01-31_CC_Spatial_SIMBA/SIMBA_in_space/data/human_DLPFC/",
    col_name=None, # column name with path to sample results
    parse_path_col=None,
    version='normal', #PCA
):
    # simba_spatial_and_RNA_SCALED_output_fn = "../results/00/simba_spatial_and_rna/mask_r500_sT_w1.0/run_simba_spatial_and_rna.output.tsv"
    results_df = pd.read_csv(results_tsv_path, sep='\t', index_col=0)

    path, fn = results_tsv_path.rsplit('/', 1)
    prefix = fn.rsplit('.', 1)[0]
    performance_fn = f"{path}/{prefix}.performance.niche.tsv"
    figures_fn = f"{path}/{prefix}.figures.niche.tsv"

    for fig_path in [f'{path}/COSSIM_CELL_TYPE_DISTR_HIST', f'{path}/COSSIM_CELL_EMB_HIST', f'{path}/COSSIM_SCATTER']:
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)

    col_name = prefix.split('.')[0] if col_name is None else col_name

    if parse_path_col is not None: # paths go directly to the adata object, rather than /path/sample/...
        results_df[col_name] = results_df[parse_path_col].apply(lambda x: x.rsplit("/", 1)[0])
        print(results_df.head())

    performance_df = pd.DataFrame(index=results_df.index)
    figure_df = pd.DataFrame(index=results_df.index)
    for s_id, r in tqdm.tqdm(results_df.iterrows(), total=results_df.shape[0]):
        sample_output_dir = r[col_name]
        adata_CG = sc.read_h5ad(f"{adata_spatial_dir}/{s_id}.h5ad")

        adata_C = sc.read_h5ad(f"{sample_output_dir}/{cell_embedding_adata_fn}")
        if version == 'PCA':
            print("setting PCA for X in new adata object")
            shared_idx = np.intersect1d(adata_C.obs.index, adata_CG.obs.index)
            adata_CG = adata_CG[shared_idx, :].copy()
            adata_C = ad.AnnData(csr_matrix(adata_C[shared_idx, :].obsm['X_pca']), obs=adata_C[shared_idx, :].obs)

        raw_data, performance_data = get_cossim_corr(adata_CG, adata_C, label='spatialLIBD', trust_n_neighbors=6)
        performance_df.loc[s_id, list(performance_data.keys())] = list(performance_data.values())

        # plot
        figure_cossim_cell_type_distr_fn = f'{path}/COSSIM_CELL_TYPE_DISTR_HIST/{s_id}.png'
        plt.hist(raw_data['cell_type_distr_cossim_flat'])
        plt.title('Cell type distribution pairwise cosine similarity')
        plt.ylabel('# pairs of cells')
        plt.xlabel('Cosine similarity of cell type distributions')
        plt.savefig(figure_cossim_cell_type_distr_fn)
        plt.close()
        figure_df.loc[s_id, 'COSSIM_CELL_TYPE_DISTR_HIST'] = figure_cossim_cell_type_distr_fn

        figure_cossim_emb_fn = f'{path}/COSSIM_CELL_EMB_HIST/{s_id}.png'
        plt.hist(raw_data['embed_cossim_flat'])
        plt.title('Embedding pairwise cosine similarity')
        plt.ylabel('# pairs of cells')
        plt.xlabel('Cosine similarity of cell embeddings')
        plt.savefig(figure_cossim_emb_fn)
        plt.close()
        figure_df.loc[s_id, 'COSSIM_CELL_EMB_HIST'] = figure_cossim_emb_fn

        figure_scatter_fn = f'{path}/COSSIM_SCATTER/{s_id}.png'
        # np.random.seed(seed=42)
        subset_idx = np.random.choice(np.arange(raw_data['embed_cossim_flat'].shape[0]), replace=False, size=10_000)
        # plot
        sns.jointplot(x=raw_data['embed_cossim_flat'][subset_idx], y=raw_data['cell_type_distr_cossim_flat'][subset_idx], linewidths=0, s=3, alpha=0.2)
        plt.xlabel('Pairwise Embedding Cosine similarity')
        plt.ylabel('Pairwise Cell type distr Cosine similarity')
        plt.savefig(figure_scatter_fn)
        plt.close()
        figure_df.loc[s_id, 'COSSIM_SCATTER'] = figure_scatter_fn

    performance_df.to_csv(performance_fn, sep='\t')
    figure_df.to_csv(figures_fn, sep='\t')
    return performance_fn, figures_fn
        