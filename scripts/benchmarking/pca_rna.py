import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import copy

import os

import argparse

def run_PCA(
    workdir, # determine which experiment is being run
    adata_CG,
):
    if not os.path.exists(workdir):
        os.m(workdir)
        
    sc.tl.pca(adata_CG)
    adata_CG.X

    adata_CG.write(os.path.join(workdir, 'adata_CG.h5ad'))
    

def main(
    workdir, # determine which experiment is being run
    adata_paths, # table of adata paths
    rerun=False
):
    output_df = pd.DataFrame(index=list(adata_paths.keys()), columns=['run_pca'])
    for sample, adata_fn in adata_paths.items():
        adata_CG = sc.read_h5ad(adata_fn)

        sample_workdir = f"{workdir}/{sample}"
        if not os.path.exists(sample_workdir) or rerun:
            run_PCA(
                workdir=sample_workdir,
                adata_CG=adata_CG,
            )

        output_df.loc[sample, 'run_pca'] = sample_workdir

    output_fn = f"{workdir}/run_pca.output.tsv"
    output_df.to_csv(output_fn, sep='\t')
    print(f"Output: {output_fn}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--workdir", type=str, help="Directory to put results")
    parser.add_argument("--adata_dir", default=None, type=str, help="Directory to anndata RNAseq object with spatial annotations")
    parser.add_argument("--adata_fns", default=None, type=str, help="list of RNAseq objects paths, comma delimited")
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
        rerun=args.rerun
    )
