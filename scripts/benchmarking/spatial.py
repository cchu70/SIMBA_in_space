import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import anndata as ad
from scipy.sparse import csr_matrix
from sklearn.metrics.pairwise import rbf_kernel
from scipy.spatial import KDTree

SPATIAL_METHODS = ['gaussian', 'rbf', 'mask']
def gen_spatial_graph(
    adata,
    e, # magnitude/mask
    scalar=1, # could be gene expression correlation?
    spatial_method='gaussian',
):
    # gaussian kernel
    spots = adata.obsm['spatial']

    if spatial_method == 'gaussian':
        squared_distances = get_squared_distances(spots)
        kernel_matrix = gaussian_kernel_matrix(squared_distances, e, scalar)
    elif spatial_method == 'rbf':
        gamma = 10**(-e)
        kernel_matrix = rbf_kernel(X=adata.obsm['spatial'], gamma=gamma)
        kernel_matrix = csr_matrix(kernel_matrix * scalar, dtype=np.float32)
    elif spatial_method == 'mask':
        kernel_matrix = get_mask(spots, radius=e)
    else:
        raise ValueError(f"spatial_method={spatial_method} not valid.")


    adata_CC = ad.AnnData(kernel_matrix)
    adata_CC.obsm['spatial'] = adata.obsm['spatial'].copy()
    adata_CC.obs = adata.obs.copy()
    adata_CC.layers['simba'] = adata_CC.X

    adata_CC.obs.index = adata.obs_names
    adata_CC.var.index = adata.obs_names
    return adata_CC

def get_mask(spots, radius=5):

    """
    ChatGPT: 
    Finds indices of neighbors within a given radius R for each point.

    :param points: List of (x, y) coordinates or a NumPy array of shape (N, 2).
    :param R: Radius within which to find neighboring points.
    :return: A list of lists, where each sublist contains indices of neighbors for a point.
    """
    tree = KDTree(spots)  # Build KDTree
    neighbors = tree.query_ball_point(spots, radius)  # Query all points for neighbors within R

    # convert to matrix
    df = pd.DataFrame(neighbors, columns=['neighbors'])
    df['value'] = 1.0
    matrix = df.explode('neighbors').pivot(columns='neighbors').fillna(0.0).to_numpy()
    return csr_matrix(matrix)

# < ChatGPT:
def get_squared_distances(spots):
    # spots: np.array list of grid points as tuples
    # Generate all grid points as (row, col) pairs
    # spots = np.array([(i, j) for i in range(N) for j in range(N)])
    # Compute the pairwise squared Euclidean distances using broadcasting
    diff = spots[:, np.newaxis, :] - spots[np.newaxis, :, :]
    squared_distances = np.sum(diff**2, axis=2)    
    return squared_distances
    
def gaussian_kernel_matrix(squared_distances, e, scalar):
    # get gamma
    max_squared_distance = np.max(squared_distances)
    gamma = max_squared_distance * 10**(-e)
    
    # Apply the Gaussian kernel
    kernel_matrix = np.exp(-squared_distances / gamma)
    # single value or nxn weight
    kernel_matrix = kernel_matrix * scalar
    return csr_matrix(kernel_matrix, dtype=np.float32)