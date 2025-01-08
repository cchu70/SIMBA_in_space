import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

from matplotlib.colors import Normalize

def plot_spatial(
    adata,
    color,  # list of genes
    x_obs_col="array_col",  # Use col as x to get the correct orientation
    y_obs_col="array_row",
    palette='viridis',
    fig_ncol=2,
    fig_size=(4,4), **kwargs
):
    # Calculate number of rows needed
    fig_nrow = int(np.ceil(len(color) / fig_ncol))
    
    # Create figure and axes
    fig, axes = plt.subplots(
        fig_nrow, fig_ncol, 
        figsize=(fig_size[0]*fig_ncol, fig_size[1]*fig_nrow), 
        sharex=True, sharey=True,
        squeeze=False  # This ensures axes is always 2D
    )
    
    # Flatten axes for easier iteration of multiple rows
    axes_flat = axes.flatten()
    
    for i, var in enumerate(color):
        ax = axes_flat[i]
        
        # Get expression values
        hues = adata[:, var].X.toarray().flatten()
        
        # Create scatter plot
        g = sns.scatterplot(
            data=adata.obs, 
            x=x_obs_col, 
            y=y_obs_col, 
            ax=ax,
            color='grey',  # We'll color the points manually
            **kwargs
        )
        
        # Create a scalar mappable for the colorbar
        # norm = Normalize(vmin=hues.min(), vmax=hues.max())
        # sm = plt.cm.ScalarMappable(cmap=palette, norm=norm)
        # sm.set_array([])
        
        # Color the points according to expression
        scatter = ax.scatter(
            adata.obs[x_obs_col],
            adata.obs[y_obs_col],
            c=hues,
            cmap=palette,
            **kwargs
        )
        
        # Add colorbar
        plt.colorbar(scatter, ax=ax)
        
        # Set title and style
        ax.set_title(var)
        ax.set_facecolor('k')
        ax.set_xlabel(None)
        ax.set_ylabel(None)
        ax.set_xticks([])
        ax.set_yticks([])
    
    # Remove empty subplots if any
    for j in range(i + 1, len(axes_flat)):
        fig.delaxes(axes_flat[j])

    ax.invert_yaxis()
    # Adjust layout to prevent overlap
    plt.tight_layout()
    
    return fig, axes
