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
from scipy.spatial import cKDTree

def combine_adatas(adata_CG, K_adata_CG, lambda_val=0.8, base_label='gene'):
    
    num_neighbors = len(K_adata_CG)
    kn_weight = np.sqrt(1 - lambda_val) / num_neighbors # uniform weight on the neighbors

    tmp_adata_CG = adata_CG.copy()
    tmp_adata_CG.X = tmp_adata_CG.X * lambda_val
    tmp_adata_CG.var['gene_type'] = base_label
    weighted_adatas = [tmp_adata_CG]

    for adata_gene_type, adata in K_adata_CG.items():
        tmp_adata = adata.copy()
        tmp_adata.X = tmp_adata.X * kn_weight
        tmp_adata.var['gene_type'] = adata_gene_type # could be rank
        weighted_adatas.append(tmp_adata)

    banksy_adata = ad.concat(weighted_adatas, axis=1)
    banksy_adata.layers['counts'] = csr_matrix(banksy_adata.X)
    return banksy_adata, lambda_val, kn_weight

def get_rankK_adata(adata_CG, k=6):

    # Build the KDTree
    coords = adata_CG.obsm['spatial']
    tree = cKDTree(coords)

    # Query the k+1 nearest neighbors (to exclude the point itself)
    distances, indices = tree.query(coords, k=k+1)
    # Exclude self (first column)
    indices = indices[:, 1:]
    
    rankK_adatas = {}
    for rank in range(k):
        k_adata = adata_CG.copy()
        k_adata.X = adata_CG.X[indices[:, rank], :]
        k_adata.var_names = k_adata.var_names + "_r" + str(rank)
        rankK_adatas[rank] = k_adata.copy()

    return rankK_adatas

method_name = 'run_simba_spatial_and_rna_rankK'

def run_simba_spatial_and_rna_rankK(
    workdir, # determine which experiment is being run
    adata_CG,
    label_col='spatialLIBD',
    k=6,
    lambda_val=0.8,
    base_label='gene',
    n_bins=20,
):
    # Set up
    si.settings.set_workdir(workdir)
    si.pp.cal_qc_rna(adata_CG)


    # Prepare graph
    rankK_adatas = get_rankK_adata(adata_CG, k=k)
    rank_adata, lambda_val, kn_weight = combine_adatas(
        adata_CG=adata_CG, 
        K_adata_CG=rankK_adatas, 
        lambda_val=lambda_val, 
        base_label=base_label
    )

    # Discretize RNA
    si.tl.discretize(rank_adata,n_bins=n_bins)
    si.pl.discretize(rank_adata,kde=False)

    si.tl.gen_graph(
        list_CG=[rank_adata],
        layer='simba',
        use_highly_variable=False, 
        dirname='graph'
    )
    # Train embedding
    si.tl.pbg_train(auto_wd=True, save_wd=True, output='model')
    
    # Read in entity embeddings obtained from pbg training.
    dict_adata = si.read_embedding()
    adata_C = dict_adata['C']  # embeddings of cells
    adata_G = dict_adata['G']
    
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
    adata_C.write(os.path.join(workdir, 'adata_C.h5ad'))
    adata_G.write(os.path.join(workdir, 'adata_G.h5ad'))

def main(
    workdir, # determine which experiment is being run
    adata_paths, # table of adata paths
    label_col='spatialLIBD',
    k=6,
    lambda_val=0.8, # weight for original gene vs neighbor mean genes
    n_bins=20,
    rerun=False,
    all=True,
    base_label='gene'
):
    output_df = pd.DataFrame(index=list(adata_paths.keys()), columns=[method_name])
    param_dir = f'k{k}_l{lambda_val}_b{n_bins}'
    
    for sample, adata_fn in adata_paths.items():
        adata_CG = sc.read_h5ad(adata_fn) 

        sample_workdir = f"{workdir}/{param_dir}/{sample}"

        adata_C_fn = f"{sample_workdir}/adata_C.h5ad"
        if not os.path.exists(adata_C_fn) or rerun:
            run_simba_spatial_and_rna_rankK(
                sample_workdir, # determine which experiment is being run
                adata_CG,
                label_col=label_col,
                k=k,
                lambda_val=lambda_val,
                base_label=base_label,
                n_bins=n_bins,
            )

        output_df.loc[sample, method_name] = sample_workdir

    if all:
        output_fn = f"{workdir}/{param_dir}/{method_name}.output.tsv"
    else:
        joined_samples = "_".join(list(adata_paths.keys()))
        output_fn = f"{workdir}/{param_dir}/{method_name}.output.{joined_samples}.tsv"
    output_df.to_csv(output_fn, sep='\t')
    print(f"Output: {output_fn}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--workdir", type=str, help="Directory to put results")
    parser.add_argument("--adata_dir", default=None, type=str, help="Directory to anndata RNAseq object with spatial annotations")
    parser.add_argument("--adata_fns", default=None, type=str, help="list of RNAseq objects paths, comma delimited")
    parser.add_argument("--label_col", default='spatialLIBD', help="Column in adata.obs table corresponding to cell labels")
    parser.add_argument("--k", default=6, type=int, help="size of nearest points for neighborhood")
    parser.add_argument("--lambda_val", default=0.8, type=float, help="value between 0 and 1 for weight of the regular gene expression")
    parser.add_argument("--n_bins", default=20, type=int, help="number of bins for discretization")
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
        workdir=args.workdir, # determine which experiment is being run
        adata_paths=adata_paths, # table of adata paths
        label_col=args.label_col,
        k=args.k,
        lambda_val=args.lambda_val, # weight for original gene vs neighbor mean genes
        n_bins=args.n_bins,
        rerun=args.rerun,
        all=args.adata_dir is not None
    )
