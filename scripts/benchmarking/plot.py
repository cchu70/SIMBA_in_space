import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import tqdm

import simba as si
import os
import scanpy as sc
import pandas as pd
import numpy as np
import copy
import seaborn as sns
from scipy.spatial import cKDTree
# import squidpy as sq

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


si.settings.set_figure_params(
    dpi=80,
    style='white',
    fig_size=[5,5],
    rc={'image.cmap': 'viridis'}
)

from matplotlib_inline.backend_inline import set_matplotlib_formats
set_matplotlib_formats('retina')

def plot_obs_spatial(
    adata, 
    obs_cols=['n_counts'],
    filter_col=None,
    filter_vals=None,
    x_obs_col="array_col",  # Use col as x to get the correct orientation
    y_obs_col="array_row",
    is_continuous=True,
    palette='viridis',
    fig_ncol=2,
    fig_size=(4,4),
    vmin=None, # list same length as obs_col
    vmax=None,
    **kwargs
):

    # Calculate number of rows needed
    fig_nrow = int(np.ceil(len(obs_cols) / fig_ncol))
    
    # Create figure and axes
    fig, axes = plt.subplots(
        fig_nrow, fig_ncol, 
        figsize=(fig_size[0]*fig_ncol, fig_size[1]*fig_nrow), 
        sharex=True, sharey=True,
        squeeze=False  # This ensures axes is always 2D
    )
    
    # Flatten axes for easier iteration of multiple rows
    axes_flat = axes.flatten()

    obs_df = copy.deepcopy(adata.obs)

    if filter_col:
        obs_df = obs_df[obs_df[filter_col].isin(filter_vals)]
    
    for i, col in enumerate(obs_cols):
        ax = axes_flat[i]

        g = sns.scatterplot(
            data=adata.obs, 
            x=x_obs_col, 
            y=y_obs_col, 
            ax=ax,
            color='lightgrey',  # We'll color the points manually
            **kwargs
        )
        
        if is_continuous:
            scatter = ax.scatter(
                x=obs_df[x_obs_col],
                y=obs_df[y_obs_col],
                c=obs_df[col],
                cmap=palette,
                vmin=obs_df[col].min() if vmin is None else vmin[i],
                vmax=obs_df[col].max() if vmax is None else vmax[i],
                **kwargs
            )
            # Add colorbar
            plt.colorbar(scatter, ax=ax)
        else:
            scatter = sns.scatterplot(
                x=obs_df[x_obs_col],
                y=obs_df[y_obs_col],
                hue=obs_df[col],
                palette=palette,
                ax=ax,
                **kwargs
            )
            plt.legend(bbox_to_anchor=(1, 1))
    
        
        # ax.set_facecolor('k')
        ax.set_xlabel(None)
        ax.set_ylabel(None)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(col)

    ax.invert_yaxis()
    return fig, ax


def plot_spatial_kernel(
    adata_CC,
    barcode_idx
):
    barcode = adata_CC.obs.index.tolist()[barcode_idx]
    fig, ax = plt.subplots()
    ax.scatter(adata_CC.obsm['spatial'][:, 0], adata_CC.obsm['spatial'][:, 1], s=3, c=adata_CC[barcode].X.toarray())
    return fig, ax

def generate_spatial_kernel_figures(
    adata_output_df,
    barcode_indices=500,# or dictionary with sample ids as the keys
    path_col="run_simba_rna_only",
    cell_embedding_adata_fn="adata_CC.h5ad",
    fig_path = "../results/00/pca_rna/PCA",
    fig_exts=["png"],
):
    fig_dir = f"{fig_path}/figures"
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)

    if isinstance(barcode_indices, int):
        barcode_indices = {s_id: barcode_indices for s_id in adata_output_df.index}
        
    for s_id, dir in tqdm.tqdm(
        adata_output_df[path_col].items(), 
        total=adata_output_df.shape[0]
    ):
        adata_CC = sc.read_h5ad(f'{dir}/{cell_embedding_adata_fn}')

        for fig_ext in fig_exts:
            fig_fn = f"{fig_path}/{s_id}.{fig_ext}"

            fig, ax = plot_spatial_kernel(adata_CC, barcode_indices[s_id])
            fig.savefig(fig_fn)
            adata_output_df.loc[s_id, f'spatial_kernel_fig_{fig_ext}'] = fig_fn

    return adata_output_df

def generate_pca_figures(
    adata_output_df,
    path_col="run_simba_rna_only",
    cell_embedding_adata_fn="adata_CG.h5ad",
    fig_path = "../results/00/pca_rna/PCA",
    adata_color_col="spatialLIBD",
    fig_exts=["png"],
    run_pca=False,
    palette=palette_celltype,
    plot_method='standard', # 'scatter_pie'
    original_adata_CG_dir=None, # path to original data for spatial coordinates. Only needed for scatter pie
    fig_size=(12, 8), 
    size=8,
):
    
    fig_dir = f"{fig_path}/figures"
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)
        
    for s_id, output_dir in tqdm.tqdm(
        adata_output_df[path_col].items(), 
        total=adata_output_df.shape[0]
    ):
        adata = sc.read_h5ad(f'{output_dir}/{cell_embedding_adata_fn}')
        if run_pca:
            si.preprocessing.pca(adata)

        for fig_ext in fig_exts:

            if plot_method == 'standard':
                fig_fn = f"{fig_path}/{s_id}.{fig_ext}"
                fig = sc.pl.pca(
                    adata, color=[adata_color_col], palette=palette, 
                    dimensions=(0, 1),
                    return_fig=True,
                    show=False,
                )
                fig.savefig(fig_fn)
                plt.close()
                adata_output_df.loc[s_id, f'pca_fig_{fig_ext}'] = fig_fn

            elif plot_method == 'scatter_pie':
                fig_fn = f"{fig_path}/{s_id}.scatter_pie.{fig_ext}"

                adata_CG = sc.read_h5ad(f'{original_adata_CG_dir}/{s_id}.h5ad')
                shared_idx = np.intersect1d(adata_CG.obs.index.to_numpy(), adata.obs.index.to_numpy())
                adata = adata[shared_idx].copy()
                adata.obsm['spatial'] = adata_CG[shared_idx, :].obsm['spatial'].copy()

                coordinates = pd.DataFrame(adata.obsm['X_pca'][:, :2], index=adata.obs.index, columns=['PC1', 'PC2'])
                plot_scatter_pie(
                    adata, 
                    coordinates, coor_x='PC1', coor_y='PC2',
                    fig_fn=fig_fn,
                    fig_size=fig_size, 
                    size=size,
                )
                plt.close()
                adata_output_df.loc[s_id, f'pca_fig_{fig_ext}'] = fig_fn


    return adata_output_df

def generate_umap_figures(
    adata_output_df,
    path_col="run_simba_rna_only",
    cell_embedding_adata_fn="adata_C.h5ad",
    original_adata_CG_dir=None, #"../data/human_DLPFC", only if plot_method == scatter_pie
    fig_path = "../results/00/simba_rna_only/UMAP",
    adata_color_col="spatialLIBD",
    fig_exts=["png"],
    run_umap=True,
    palette=palette_celltype,
    include_legend=True,
    plot_method='standard', # ['scatter_pie']
    fig_size=(12, 8), 
    size=8,
):
    """
    adata_output_df: 
        index: sample names
    path_col: 
        column in adata_output_df with paths to the output e.g. "./results/00/simba_rna_only/151675"
    fig_path:
        path to save figures
    adata_color_col:
        column in adata.obs that 
    fig_exts:
        png, svg, etc.
    """
    si.settings.set_workdir(fig_path)
    
    for s_id, dir in tqdm.tqdm(
        adata_output_df[path_col].items(), 
        total=adata_output_df.shape[0]
    ):  

        adata_C = sc.read_h5ad(f'{dir}/{cell_embedding_adata_fn}')

        if run_umap:
            si.tl.umap(adata_C,n_neighbors=15,n_components=2)

        for fig_ext in fig_exts:

            if plot_method == 'standard':
                fig_fn = f"{s_id}.{fig_ext}"
                si.pl.umap(
                    adata_C,color=[adata_color_col],
                    dict_palette={adata_color_col: palette} if palette is not None else None,
                    fig_size=(6,4),
                    drawing_order='random',
                    save_fig=True,
                    fig_name=fig_fn,
                    legend="auto" if include_legend else False,
                )
                adata_output_df.loc[s_id, f'umap_fig_{fig_ext}'] = f"{fig_path}/figures/{fig_fn}"
            elif plot_method == 'scatter_pie':
                fig_fn = f"{s_id}.scatter_pie.{fig_ext}"

                adata_CG = sc.read_h5ad(f'{original_adata_CG_dir}/{s_id}.h5ad')
                adata_C.obsm['spatial'] = adata_CG[adata_C.obs.index, :].obsm['spatial'].copy()

                coordinates = pd.DataFrame(adata_C.obsm['X_umap'], index=adata_C.obs.index, columns=['umap1', 'umap2'])
                plot_scatter_pie(
                    adata_C, 
                    coordinates, coor_x='umap1', coor_y='umap2',
                    fig_size=fig_size,
                    size=size,
                    fig_fn=f'{fig_path}/{fig_fn}'
                )
                adata_output_df.loc[s_id, f'umap_fig_{fig_ext}'] = f'{fig_path}/{fig_fn}'

            plt.close()
            

    return adata_output_df

def draw_pie(dist, 
             xpos, 
             ypos, 
             size, 
             colors,
             ax=None):
    """
    From here: https://stackoverflow.com/a/56338575
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10,8))

    # for incremental pie slices
    non_zero_dist = dist[dist > 0]
    cumsum = np.cumsum(non_zero_dist)
    cumsum = cumsum/ cumsum[-1]
    pie = [0] + cumsum.tolist()

    i = 0
    for r1, r2 in zip(pie[:-1], pie[1:]):
        angles = np.linspace(2 * np.pi * r1, 2 * np.pi * r2)
        x = [0] + np.cos(angles).tolist()
        y = [0] + np.sin(angles).tolist()

        xy = np.column_stack([x, y])

        label = non_zero_dist.index[i]
        ax.scatter([xpos], [ypos], marker=xy, s=size, facecolor=colors[label], linewidths=0.1, edgecolors='k')
        i += 1
    
    return ax

def plot_scatter_pie(adata_C, coordinates, coor_x='umap1', coor_y='umap2', fig_size=(4, 4), size=8, fig_fn=None):

    coords = adata_C.obsm['spatial']
    tree = cKDTree(coords)

    k=6
    # Query the k+1 nearest neighbors (to exclude the point itself)
    distances, indices = tree.query(coords, k=k+1)
    # Exclude self (first column)
    indices = indices[:, 1:]
    indices_df = pd.DataFrame(indices, index=adata_C.obs.index, columns=np.arange(1, k + 1))


    fig, ax = plt.subplots(figsize=fig_size)
    for cell, r in coordinates.iterrows():
        k_neighbors = indices_df.loc[cell].to_numpy()
        dist = adata_C.obs.loc[adata_C.obs.index[k_neighbors]]['spatialLIBD'].value_counts()
        draw_pie(dist, 
                r[coor_x], 
                r[coor_y], 
                size=size, 
                colors=palette_celltype,
                ax=ax)
        
    fig.savefig(fig_fn)
    plt.close()

def generate_spatial_figures(
    adata_output_df,
    path_col="run_simba_rna_only",
    cell_embedding_adata_fn="adata_C.h5ad",
    adata_spatial_dir="/data/pinello/PROJECTS/2025-01-31_CC_Spatial_SIMBA/SIMBA_in_space/data/human_DLPFC/", # directory with the spatial coordinates.e.g.  <sample_id>.h5ad
    fig_path = "../results/00/simba_rna_only/SPATIAL",
    adata_color_col="spatialLIBD", # leiden, etc.
    fig_exts=["png"],
    **kwargs,
):
    if not os.path.exists(fig_path):
        os.makedirs(fig_path)

    for s_id, dir in tqdm.tqdm(
        adata_output_df[path_col].items(), 
        total=adata_output_df.shape[0]
    ):
        orig_adata = sc.read_h5ad(f"{adata_spatial_dir}/{s_id}.h5ad")
        adata = sc.read_h5ad(f'{dir}/{cell_embedding_adata_fn}')
        adata.obsm['spatial'] = orig_adata[adata.obs.index].obsm['spatial'].copy()
        adata.obs['array_col'] = orig_adata[adata.obs.index].obsm['spatial'].copy()[:, 0]
        adata.obs['array_row'] = orig_adata[adata.obs.index].obsm['spatial'].copy()[:, 1]

        for fig_ext in fig_exts:
            fig_fn = f"{fig_path}/{s_id}.{fig_ext}"
            fig, ax = plot_obs_spatial(
                adata,
                obs_cols=[adata_color_col],
                palette=None,
                is_continuous=False,
                **kwargs
            )
            fig.savefig(fig_fn)
            plt.close()
            adata_output_df.loc[s_id, f'spatial_fig_{fig_ext}'] = fig_fn

    return adata_output_df
    
    

def combine_images(
    adata_output_df,
    img_path_col="umap_fig_png",
    nrows=3,
    ncols=4,
    title=None,
    x_min=110, x_max=950,
    y_min=100, y_max=740
):
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, 10), sharex=True, sharey=True)
    axes_flat = axes.flatten()
    i = 0
    for s_id, fig_fn in adata_output_df[img_path_col].items():
        img = mpimg.imread(fig_fn)
        axes_flat[i].imshow(img[y_min:y_max, x_min:x_max, :], interpolation='nearest', aspect='equal') # (840, 1185, 4)
        axes_flat[i].set_title(s_id)
        axes_flat[i].set_axis_off()
        i += 1
    fig.suptitle(title)
    plt.xticks([])
    plt.yticks([])
    plt.tight_layout()
    return fig, axes
    
        