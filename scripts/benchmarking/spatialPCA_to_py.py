import anndata as ad
import pandas as pd
import numpy as np
import os
import argparse
import scanpy as sc

def convert_to_adata(spatialpca_sample_dir):
    spatialPCs_fn = f"{spatialpca_sample_dir}/Spatial_PCA.spatialPCs.tsv"
    counts_fn = f"{spatialpca_sample_dir}/Spatial_PCA.counts.mtx"
    features_fn = f"{spatialpca_sample_dir}/Spatial_PCA.counts_features.tsv"
    # barcodes_fn = f"{spatialpca_sample_dir}/Spatial_PCA.counts_barcodes.tsv"
    obs_fn = f"{spatialpca_sample_dir}/Spatial_PCA.obs_table.tsv"
    # norm_expr_fn = f"{spatialpca_sample_dir}/Spatial_PCA.normalized_expr.mtx"

    spatialPCs_df = pd.read_csv(spatialPCs_fn, sep='\t').T # saved as d x barcodes
    

    adata = sc.read_mtx(counts_fn).T 
    features = pd.read_csv(features_fn, sep='\t', index_col=0, header=None).index
    obs = pd.read_csv(obs_fn, sep='\t')
    obs = obs.set_index('barcode') # in some samples, there is a 1.1 suffix added in the index
    adata.obs = obs.copy() #.index = barcodes
    adata.var.index = features
    adata.var.index.name = None

    # Not the same shape as counts?
    # norm_expr_adata = sc.read_mtx(norm_expr_fn).T 
    # norm_expr_adata.obs.index = barcodes
    # norm_expr_adata.var.index = features
    # adata.layers['norm'] = norm_expr_adata.X

    # some runs had some barcodes filtered out before getting the PCs
    PC_adata = adata[spatialPCs_df.index, :].copy()

    PC_adata.obsm['X_pca'] = spatialPCs_df.to_numpy()
    output_fn = os.path.join(spatialpca_sample_dir, 'spatialPCA_to_py.adata.h5ad')
    print(output_fn)
    PC_adata.write(output_fn)

    return output_fn

def main(spatialpca_dir, samples):
    spatialpca_to_py_fn = f"{spatialpca_dir}/spatialpca_adata_py.tsv"
    spatialpca_to_py_df = pd.DataFrame()

    samples = samples.split(",") if samples is not None else os.listdir(spatialpca_dir)
    for sample in samples:
        full_sample_dir = f"{spatialpca_dir}/{sample}"
        if os.path.isdir(full_sample_dir):
            adata_fn = convert_to_adata(full_sample_dir)
            spatialpca_to_py_df.loc[sample, 'spatialPCA_adata_fn'] = adata_fn

    spatialpca_to_py_df.to_csv(spatialpca_to_py_fn, sep='\t')
    
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--spatialpca_dir", default=None, type=str, help="Directory each sample's results from SpatialPCA (spatialpca/sample1, etc..)")
    parser.add_argument("--samples", default=None, type=str, help="comma delimited list of subdirectories in --spatialpca_dir corresponding to sample names")
    
    args = parser.parse_args()
    main(args.spatialpca_dir, args.samples)