import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import copy

import os
import simba as si

import argparse

def run_simba_rna_only(
    workdir, # determine which experiment is being run
    adata_CG,
    label_col='spatialLIBD'
):
    # Set up
    si.settings.set_workdir(workdir)
    si.pp.cal_qc_rna(adata_CG)
    
    # Prepare graph
    si.tl.discretize(adata_CG,n_bins=5)
    si.pl.discretize(adata_CG,kde=False)

    si.tl.gen_graph(
        list_CG=[adata_CG],
        layer='simba',
        use_highly_variable=False, 
        dirname='graph0'
    )

    # Train embedding
    si.tl.pbg_train(auto_wd=True, save_wd=True, output='model')
    
    # Read in entity embeddings obtained from pbg training.
    dict_adata = si.read_embedding()
    adata_C = dict_adata['C']  # embeddings of cells
    adata_G = dict_adata['G']  # embeddings of genes
    
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
    adata_CG.write(os.path.join(workdir, 'adata_CG.h5ad'))
    adata_C.write(os.path.join(workdir, 'adata_C.h5ad'))
    adata_G.write(os.path.join(workdir, 'adata_G.h5ad'))

def main(
    workdir, # determine which experiment is being run
    adata_paths, # table of adata paths
    label_col='spatialLIBD',
    rerun=False
):
    output_df = pd.DataFrame(index=list(adata_paths.keys()), columns=['run_simba_rna_only'])
    for sample, adata_fn in adata_paths.items():
        adata_CG = sc.read_h5ad(adata_fn)

        sample_workdir = f"{workdir}/{sample}"
        if not os.path.exists(sample_workdir) or rerun:
            run_simba_rna_only(
                workdir=sample_workdir,
                adata_CG=adata_CG,
                label_col=label_col,
            )

        output_df.loc[sample, 'run_simba_rna_only'] = sample_workdir

    output_fn = f"{workdir}/run_simba_rna_only.output.tsv"
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
