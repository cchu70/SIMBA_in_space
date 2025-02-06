import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import tqdm

import simba as si
import os
import scanpy as sc
import pandas as pd

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

def generate_pca_figures(
    adata_output_df,
    path_col="run_simba_rna_only",
    cell_embedding_adata_fn="adata_CG.h5ad",
    fig_path = "../results/00/pca_rna/PCA",
    adata_color_col="spatialLIBD",
    fig_exts=["png"],
    run_pca=False,
):
    
    fig_dir = f"{fig_path}/figures"
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)
        
    for s_id, dir in tqdm.tqdm(
        adata_output_df[path_col].items(), 
        total=adata_output_df.shape[0]
    ):
        adata = sc.read_h5ad(f'{dir}/{cell_embedding_adata_fn}')
        
        if run_pca:
            si.preprocessing.pca(adata)

        for fig_ext in fig_exts:
            fig_fn = f"{fig_path}/{s_id}.{fig_ext}"
            fig = sc.pl.pca(
                adata, color=[adata_color_col], palette=palette_celltype, 
                dimensions=(0, 1),
                return_fig=True,
                show=False,
            )
            fig.savefig(fig_fn)
            adata_output_df.loc[s_id, f'pca_fig_{fig_ext}'] = fig_fn

    return adata_output_df

def generate_umap_figures(
    adata_output_df,
    path_col="run_simba_rna_only",
    cell_embedding_adata_fn="adata_C.h5ad",
    fig_path = "../results/00/simba_rna_only/UMAP",
    adata_color_col="spatialLIBD",
    fig_exts=["png"],
    run_umap=True,
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
            fig_fn = f"{s_id}.{fig_ext}"
            si.pl.umap(
                adata_C,color=[adata_color_col],
                dict_palette={adata_color_col: palette_celltype},
                fig_size=(6,4),
                drawing_order='random',
                save_fig=True,
                fig_name=fig_fn
            )
            adata_output_df.loc[s_id, f'umap_fig_{fig_ext}'] = f"{fig_path}/figures/{fig_fn}"

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
    
        