import os
import pandas as pd
import numpy as np
import json
from anndata import (
    AnnData,
    read_h5ad,
    read_csv,
    read_excel,
    read_hdf,
    read_loom,
    read_mtx,
    read_text,
    read_umi_tools,
    read_zarr,
)
from pathlib import Path
import tables

def download_spatial_dlpfc(
    sample_id,
    workdir,
    file_type="h5_filtered",
    link_prefix="https://spatial-dlpfc.s3.us-east-2.amazonaws.com/",
):
    """
    sample_id: 
        sample id. See valid IDs at https://research.libd.org/spatialLIBD/
    file_type:
        h5_filtered or h5_raw
    """

    file_type_suffixes = {
        "h5_raw": f"h5/{sample_id}_raw_feature_bc_matrix.h5",
        "h5_filtered": f"h5/{sample_id}_filtered_feature_bc_matrix.h5"
    }

    download_path = f"{link_prefix}/{file_type_suffixes[file_type]}"
    output_path = f"{workdir}/data"

    if not os.path.exists(output_path):
        os.makedirs(output_path)
    
    print(f"download_path={download_path} to {output_path}")
    os.system(f'wget "{download_path}"')
    
    return download_path

def read_v3_10x_h5(filename):
    """
    Read hdf5 file from Cell Ranger v3 or later versions.
    """
    with tables.open_file(str(filename), 'r') as f:
        try:
            dsets = {}
            for node in f.walk_nodes('/', 'Array'): # for node in f.walk_nodes('/matrix', 'Array'):
                dsets[node.name] = node.read()
            from scipy.sparse import csr_matrix

            print(dsets.keys())
            print(dsets['data'].shape)
            print(dsets['indices'])
            print(dsets['indptr'])

            # M, N = dsets['shape']
            data = dsets['data']
            if dsets['data'].dtype == np.dtype('int32'):
                data = dsets['data'].view('float32')
                data[:] = dsets['data']
            matrix = csr_matrix(
                (data, dsets['indices'], dsets['indptr']),
                # shape=(N, M),
            )
            adata = AnnData(
                matrix,
                obs=dict(obs_names=dsets['barcodes'].astype(str)),
                var=dict(
                    var_names=dsets['name'].astype(str),
                    gene_ids=dsets['id'].astype(str),
                    feature_types=dsets['feature_type'].astype(str),
                    genome=dsets['genome'].astype(str),
                ),
            )
            return adata
        except KeyError:
            raise Exception('File is missing one or more required datasets.')

def read_legacy_10x_h5(filename, genome=None):
    """
    Read hdf5 file from Cell Ranger v2 or earlier versions.
    """
    with tables.open_file(str(filename), 'r') as f:
        try:
            children = [x._v_name for x in f.list_nodes(f.root)]
            if not genome:
                if len(children) > 1:
                    raise ValueError(
                        f"'{filename}' contains more than one genome. "
                        "For legacy 10x h5 "
                        "files you must specify the genome "
                        "if more than one is present. "
                        f"Available genomes are: {children}"
                    )
                genome = children[0]
            elif genome not in children:
                raise ValueError(
                    f"Could not find genome '{genome}' in '{filename}'. "
                    f'Available genomes are: {children}'
                )
            dsets = {}
            for node in f.walk_nodes('/' + genome, 'Array'):
                dsets[node.name] = node.read()
            # AnnData works with csr matrices
            # 10x stores the transposed data, so we do the transposition
            from scipy.sparse import csr_matrix

            M, N = dsets['shape']
            data = dsets['data']
            if dsets['data'].dtype == np.dtype('int32'):
                data = dsets['data'].view('float32')
                data[:] = dsets['data']
            matrix = csr_matrix(
                (data, dsets['indices'], dsets['indptr']),
                shape=(N, M),
            )
            # the csc matrix is automatically the transposed csr matrix
            # as scanpy expects it, so, no need for a further transpostion
            adata = AnnData(
                matrix,
                obs=dict(obs_names=dsets['barcodes'].astype(str)),
                var=dict(
                    var_names=dsets['gene_names'].astype(str),
                    gene_ids=dsets['genes'].astype(str),
                ),
            )
            return adata
        except KeyError:
            raise Exception('File is missing one or more required datasets.')
            
