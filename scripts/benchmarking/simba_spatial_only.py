import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import copy

import os
import simba as si

import argparse

def gen_spatial_graph(
    adata,
    x_col,
    y_col,
    gamma,
):
    # gaussian kernel
    spots = adata.obs[[x_col, y_col]].toarray()

    kernel_matrix = gaussian_kernel_matrix(spots, gamma)
    pass

# ChatGPT:
def gaussian_kernel_matrix(spots, gamma):

    # spots: np.array list of grid points as tuples
    # Generate all grid points as (row, col) pairs
    # spots = np.array([(i, j) for i in range(N) for j in range(N)])
    
    # Compute the pairwise squared Euclidean distances using broadcasting
    diff = spots[:, np.newaxis, :] - spots[np.newaxis, :, :]
    squared_distances = np.sum(diff**2, axis=2)
    
    # Apply the Gaussian kernel
    kernel_matrix = np.exp(-squared_distances / gamma)
    
    return kernel_matrix

def run_simba_spatial_only(
    workdir, # determine which experiment is being run
    adata_CG,
    label_col='spatialLIBD'
):
    # Set up
    si.settings.set_workdir(workdir)
    si.pp.cal_qc_rna(adata_CG)
    
    # Prepare graph
    adata_CC = None

    si.tl.gen_graph(
        list_adata=[adata_CC],
        prefix='C',
        layer='simba',
        use_highly_variable=False, 
        dirname='graph0'
    )

    # Train embedding
    si.tl.pbg_train(auto_wd=True, save_wd=True, output='model')
    
    # Read in entity embeddings obtained from pbg training.
    dict_adata = si.read_embedding()
    adata_C = dict_adata['C']  # embeddings of cells
    
    adata_C.obs[label_col] = adata_CG[adata_C.obs_names,:].obs[label_col].copy()
    adata_C.obs['n_counts'] = adata_CG[adata_C.obs_names,:].obs['n_counts'].copy()
    adata_C.obs['n_genes'] = adata_CG[adata_C.obs_names,:].obs['n_genes'].copy()
    
    # save
    adata_CC.write(os.path.join(workdir, 'adata_CC.h5ad'))
    adata_C.write(os.path.join(workdir, 'adata_C.h5ad'))

def main(
    workdir, # determine which experiment is being run
    adata_paths, # table of adata paths
    label_col='spatialLIBD',
    rerun=False
):
    output_df = pd.DataFrame(index=list(adata_paths.keys()), columns=['run_simba_spatial_only'])
    for sample, adata_fn in adata_paths.items():
        adata_CG = sc.read_h5ad(adata_fn)

        sample_workdir = f"{workdir}/{sample}"
        if not os.path.exists(sample_workdir) or rerun:
            run_simba_spatial_only(
                workdir=sample_workdir,
                adata_CG=adata_CG,
                label_col=label_col,
            )

        output_df.loc[sample, 'run_simba_spatial_only'] = sample_workdir

    output_fn = f"{workdir}/run_simba_spatial_only.output.tsv"
    output_df.to_csv(output_fn, sep='\t')
    print(f"Output: {output_fn}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--workdir", type=str, help="Directory to put results")
    parser.add_argument("--adata_dir", default=None, type=str, help="Directory to anndata RNAseq object with spatial annotations")
    parser.add_argument("--adata_fns", default=None, type=str, help="list of RNAseq objects paths, comma delimited")
    parser.add_argument("--label_col", default='spatialLIBD', help="Column in adata.obs table corresponding to cell labels")
    parser.add_argument("--rerun", action=argparse.BooleanOptionalAction, default=False, help="Rerun")
    
    args = parser.parse_args()

    adata_paths = {}

    if args.adata_dir:
        fn_list = os.listdir(args.adata_dir)
    elif args.adata_fns:
        fn_list = args.adata_fns.split(",")
         
    for fn in os.listdir(args.adata_dir):
        sample = fn.rsplit('.', 1)[0]
        adata_paths[sample] = f"{args.adata_dir}/{fn}"

    if not os.path.exists(args.workdir):
        print("Making directory")
        os.makedirs(args.workdir)

    print(f"Running {len(adata_paths)} files.")
    main(
        workdir=args.workdir,
        adata_paths=adata_paths,
        label_col=args.label_col,
        rerun=args.rerun
    )
