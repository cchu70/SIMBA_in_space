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
from spatial import gen_spatial_graph, SPATIAL_METHODS

def run_simba_spatial_only(
    workdir, # determine which experiment is being run
    adata_CG,
    e = 2, # fraction of max distance for gamma
    spatial_method='gaussian',
    label_col='spatialLIBD'
):
    # Set up
    si.settings.set_workdir(workdir)
    si.pp.cal_qc_rna(adata_CG)
    
    # Prepare graph
    adata_CC = gen_spatial_graph(adata_CG, e = e, spatial_method=spatial_method)

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
    adata_C = dict_adata['C0']  # embeddings of cells
    
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
    spatial_method='gaussian',
    e=2,
    rerun=False
):
    output_df = pd.DataFrame(index=list(adata_paths.keys()), columns=['run_simba_spatial_only'])

    assert spatial_method in SPATIAL_METHODS

    if spatial_method in ['gaussian', 'rbf']:
        spatial_method_dir = f'{spatial_method}_e{e}' 
    elif spatial_method == 'mask':
        spatial_method_dir = f"{spatial_method}_r{e}" # mask radius

    for sample, adata_fn in adata_paths.items():
        adata_CG = sc.read_h5ad(adata_fn)

        sample_workdir = f"{workdir}/{spatial_method_dir}/{sample}"
        print(f"{sample_workdir}")
        if not os.path.exists(sample_workdir) or rerun:
            run_simba_spatial_only(
                workdir=sample_workdir,
                adata_CG=adata_CG,
                label_col=label_col,
                spatial_method=spatial_method,
                e=e,
            )

        output_df.loc[sample, 'run_simba_spatial_only'] = sample_workdir

    output_fn = f"{workdir}/{spatial_method_dir}/run_simba_spatial_only.output.tsv"
    output_df.to_csv(output_fn, sep='\t')
    print(f"Output: {output_fn}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--workdir", type=str, help="Directory to put results")
    parser.add_argument("--adata_dir", default=None, type=str, help="Directory to anndata RNAseq object with spatial annotations")
    parser.add_argument("--adata_fns", default=None, type=str, help="list of RNAseq objects paths, comma delimited")
    parser.add_argument("--label_col", default='spatialLIBD', help="Column in adata.obs table corresponding to cell labels")
    parser.add_argument("--spatial_method", default='gaussian', help='Approach for spatial edges: gaussian, rbf, or mask')
    parser.add_argument("--e", default=2, type=int, help="Parameter for corresponding spatial method (gamma for gaussian/rbf, radius for mask)")
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
        e=args.e,
        spatial_method=args.spatial_method,
        rerun=args.rerun
    )
