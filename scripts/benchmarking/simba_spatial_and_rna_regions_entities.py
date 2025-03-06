import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import copy

import os
import simba as si

import argparse

import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import anndata as ad
from scipy.sparse import csr_matrix
from sklearn.metrics.pairwise import rbf_kernel
from spatial import gen_spatial_graph, SPATIAL_METHODS
import squidpy as sq
# ChatGPT>    

def gen_grid(adata, grid_size):
    max_x = np.max(adata.obsm['spatial'][:, 0])
    max_y = np.max(adata.obsm['spatial'][:, 1])

    min_x = np.min(adata.obsm['spatial'][:, 0])
    min_y = np.min(adata.obsm['spatial'][:, 1])

    x_grid = np.arange(max(0, min_x - grid_size), max_x + grid_size, step=grid_size)
    y_grid = np.arange(max(0, min_y - grid_size), max_y + grid_size, step=grid_size)

    for i in range(len(x_grid) - 1):
        col_fil = (adata.obsm['spatial'][:, 0] >= x_grid[i]) & (adata.obsm['spatial'][:, 0] < x_grid[i + 1])
        adata.obs.loc[col_fil, 'x_region'] = i
        adata.obs.loc[col_fil, 'x_region_val'] = x_grid[i]
        
    for j in range(len(y_grid) - 1):
        row_fil = (adata.obsm['spatial'][:, 1] >= y_grid[j]) & (adata.obsm['spatial'][:, 1] < y_grid[j + 1])
        adata.obs.loc[row_fil, 'y_region'] = j
        adata.obs.loc[row_fil, 'y_region_val'] = y_grid[j]

    adata.obs['region_id'] = adata.obs['x_region'].astype(int).astype(str) + "_" + adata.obs['y_region'].astype(int).astype(str)
    adata.obs['region_val_id'] = adata.obs['x_region_val'].astype(int).astype(str) + "_" + adata.obs['y_region_val'].astype(int).astype(str)
    adata.uns['grid'] = {}
    adata.uns['grid']['x_grid'] = x_grid
    adata.uns['grid']['y_grid'] = y_grid
    return adata

def gen_cell_region(adata):
    """
    Generate cell x region adata object
    """
    cell_region_obs_df = adata.obs[
        ['x_region', 'y_region', 'region_id']
        ].reset_index().rename(columns={'index': 'barcode'}).sort_values(by=['x_region', 'y_region'])
    cell_region_obs_df['value'] = 1 # 1 if the cell is in the region
    cell_region_X = cell_region_obs_df.pivot(
        index='barcode', columns='region_id', values='value'
    ).fillna(0)

    cell_region_var_df = cell_region_X.sum().to_frame().rename(columns={0: "num_barcodes"})

    adata_CR = ad.AnnData(csr_matrix(cell_region_X.to_numpy()), obs=adata.obs, var=cell_region_var_df)
    adata_CR.layers['simba'] = adata_CR.X.copy()
    return adata_CR
            

def plot_spatial_grid(adata, color='MOBP', scale_factor=1):
    """
    scale_factor: adata.uns['spatial']['151507']['scalefactors']['tissue_hires_scalef']
    """
    ax = sq.pl.spatial_scatter(adata, color=color, size=1, linewidth=0, return_ax=True)
    x_grid = adata.uns['grid']['x_grid']
    y_grid = adata.uns['grid']['y_grid']
    for x in x_grid:
        ax.axvline(x * scale_factor, c='k') 

    for y in y_grid:
        ax.axhline(y * scale_factor, c='k') 

    return ax

def run_simba_spatial_and_rna_region_entitites(
    workdir, # determine which experiment is being run
    adata_CG,
    label_col='spatialLIBD',
    grid_size=500,
    CR_training_weight=1.0 # weight for CR edge type.
):
    # Set up
    si.settings.set_workdir(workdir)
    si.pp.cal_qc_rna(adata_CG)


    # Prepare graph
    adata_CG = gen_grid(adata_CG, grid_size)
    adata_CR = gen_cell_region(adata_CG)

    # Discretize RNA
    si.tl.discretize(adata_CG,n_bins=5)
    si.pl.discretize(adata_CG,kde=False)

    si.tl.gen_graph(
        list_adata=[adata_CG, adata_CR],
        add_edge_weights=True,
        layer='simba',
        use_highly_variable=False, 
        dirname='graph'
    )

    # Update weight
    reweight_CR_dict_config = si.settings.pbg_params.copy()
    # Get the C-C edge type param
    CR_weight_idx = [i for i, param in enumerate(reweight_CR_dict_config['relations']) if (param['lhs'] == 'E0') and (param['rhs'] == 'E2')]
    assert len(CR_weight_idx) == 1
    
    reweight_CR_dict_config['relations'][CR_weight_idx[0]]['weight'] = CR_training_weight

    # Train embedding
    si.tl.pbg_train(pbg_params = reweight_CR_dict_config, auto_wd=True, save_wd=True, output='model')
    
    # Read in entity embeddings obtained from pbg training.
    dict_adata = si.read_embedding()
    adata_C = dict_adata['E0']  
    adata_G = dict_adata['E1']  
    adata_R = dict_adata['E2']
    
    adata_C.obs[label_col] = adata_CG[adata_C.obs_names,:].obs[label_col].copy()
    adata_C.obs['n_counts'] = adata_CG[adata_C.obs_names,:].obs['n_counts'].copy()
    adata_C.obs['n_genes'] = adata_CG[adata_C.obs_names,:].obs['n_genes'].copy()
    
    # embed cells and genes into the same space
    adata_all = si.tl.embed(adata_ref=adata_C,list_adata_query=[adata_G])

    # Compare entities
    adata_cmp = si.tl.compare_entities(adata_ref=adata_C, adata_query=adata_G)
    
    # save
    adata_all.write(os.path.join(workdir, 'adata_all.h5ad'))
    adata_cmp.write(os.path.join(workdir, 'adata_cmp.h5ad'))
    adata_R.write(os.path.join(workdir, 'adata_R.h5ad'))
    adata_C.write(os.path.join(workdir, 'adata_C.h5ad'))
    adata_G.write(os.path.join(workdir, 'adata_G.h5ad'))
    adata_CG.write(os.path.join(workdir, 'adata_CG.h5ad'))
    adata_CR.write(os.path.join(workdir, 'adata_CR.h5ad'))

def main(
    workdir, # determine which experiment is being run
    adata_paths, # table of adata paths
    label_col='spatialLIBD',
    grid_size=500,
    CR_training_weight=1.0,
    rerun=False,
    all=True
):
    output_df = pd.DataFrame(index=list(adata_paths.keys()), columns=['run_simba_spatial_and_rna']) # typo

    param_dir = f'regions_g{grid_size}_w{CR_training_weight}' 

    for sample, adata_fn in adata_paths.items():
        adata_CG = sc.read_h5ad(adata_fn)

        sample_workdir = f"{workdir}/{param_dir}/{sample}"
        if not os.path.exists(sample_workdir) or rerun:
            run_simba_spatial_and_rna_region_entitites(
                workdir=sample_workdir,
                adata_CG=adata_CG,
                label_col=label_col,
                grid_size=grid_size,
                CR_training_weight=CR_training_weight,
            )

        output_df.loc[sample, 'run_simba_spatial_and_rna'] = sample_workdir

    if all:
        output_fn = f"{workdir}/run_simba_spatial_and_rna.output.tsv"
    else:
        joined_samples = "_".join(list(adata_paths.keys()))
        output_fn = f"{workdir}/run_simba_spatial_and_rna.output.{joined_samples}.tsv"
    output_df.to_csv(output_fn, sep='\t')
    print(f"Output: {output_fn}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--workdir", type=str, help="Directory to put results")
    parser.add_argument("--adata_dir", default=None, type=str, help="Directory to anndata RNAseq object with spatial annotations")
    parser.add_argument("--adata_fns", default=None, type=str, help="list of RNAseq objects paths, comma delimited")
    parser.add_argument("--label_col", default='spatialLIBD', help="Column in adata.obs table corresponding to cell labels")
    parser.add_argument("--grid_size", default=500, type=float, help="Size of square grids to define regions on the slide")
    parser.add_argument("--CR_training_weight", default=1, type=float, help="cell-region edge type weight during training.")
    parser.add_argument("--rerun", action=argparse.BooleanOptionalAction, default=False, help="Rerun")
    
    args = parser.parse_args()

    adata_paths = {}

    if args.adata_dir:
        fn_list = os.listdir(args.adata_dir) # returns file name only
        for fn in fn_list:
            sample = fn.rsplit('.', 1)[0]
            adata_paths[sample] = f"{args.adata_dir}/{fn}"
            
    elif args.adata_fns:
        fn_list = args.adata_fns.split(",")
        for fn in fn_list:
            sample = fn.split('/')[-1].rsplit('.', 1)[0]
            adata_paths[sample] = fn

    if not os.path.exists(args.workdir):
        print("Making directory")
        os.makedirs(args.workdir)

    print(f"Running {len(adata_paths)} files.")
    main(
        workdir=args.workdir,
        adata_paths=adata_paths,
        label_col=args.label_col,
        grid_size=args.grid_size,
        CR_training_weight=args.CR_training_weight,
        rerun=args.rerun,
        all=args.adata_dir is not None
    )
