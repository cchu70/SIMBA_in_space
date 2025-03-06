import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import copy
from matplotlib.colors import Normalize

# approximate original figure in http://spatial.libd.org/spatialLIBD/
palette_celltype={'L1':'#eb34a8',
                  'L2':'#3486eb',
                  'L3':'#34eb5b',
                  'L4':"#ae34eb",
                  'L5':'#ebdb34',
                  'L6':'#eb9234',
                  'WM':'#000000'}

palette_entity_anno = palette_celltype.copy()
palette_entity_anno['gene'] = "lightgray"

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


def plot_obs_spatial(
    adata, 
    obs_cols=['n_counts'],
    filter_col=None,
    filter_vals=None,
    x_obs_col="array_col",  # Use col as x to get the correct orientation
    y_obs_col="array_row",
    palette='viridis',
    fig_ncol=2,
    fig_size=(4,4),
    vmin=None, # list same length as obs_col
    vmax=None,
    **kwargs
):

    # Calculate number of rows needed
    fig_nrow = int(np.ceil(len(obs_cols) / fig_ncol))
    
    # Create figure and axes
    fig, axes = plt.subplots(
        fig_nrow, fig_ncol, 
        figsize=(fig_size[0]*fig_ncol, fig_size[1]*fig_nrow), 
        sharex=True, sharey=True,
        squeeze=False  # This ensures axes is always 2D
    )
    
    # Flatten axes for easier iteration of multiple rows
    axes_flat = axes.flatten()

    obs_df = copy.deepcopy(adata.obs)

    if filter_col:
        obs_df = obs_df[obs_df[filter_col].isin(filter_vals)]
    
    for i, col in enumerate(obs_cols):
        ax = axes_flat[i]

        g = sns.scatterplot(
            data=adata.obs, 
            x=x_obs_col, 
            y=y_obs_col, 
            ax=ax,
            color='lightgrey',  # We'll color the points manually
            **kwargs
        )
        
        scatter = ax.scatter(
            x=obs_df[x_obs_col],
            y=obs_df[y_obs_col],
            c=obs_df[col],
            cmap=palette,
            vmin=obs_df[col].min() if vmin is None else vmin[i],
            vmax=obs_df[col].max() if vmax is None else vmax[i],
            **kwargs
        )
    
        # Add colorbar
        plt.colorbar(scatter, ax=ax)
        
        ax.set_facecolor('k')
        ax.set_xlabel(None)
        ax.set_ylabel(None)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(col)

    ax.invert_yaxis()
    return fig, ax

# from Yo Akiyama
import matplotlib.ticker as ticker
import scipy.stats as stats
def setup_figure(aw=4.5, ah=3, xspace=[0.75,0.25], yspace=[0.75,0.25],
                 colorbar=False, ds=0.15, cw=0.15, ct=0, ch=None):
    dl, dr = xspace
    db, dt = yspace
    fw = dl + aw + dr
    fh = db + ah + dt
    fig = plt.figure(facecolor=(1,1,1), figsize=(fw,fh))
    ax = fig.add_axes([dl/fw, db/fh, aw/fw, ah/fh])
    if not colorbar:
        return ax
    else:
        if ch is None:
            ch = ah/2
        cax = fig.add_axes([(dl+aw+ds)/fw, (db+ah-ch-ct)/fh, cw/fw, ch/fh])
        return ax, cax

def qqplot(pval, pval_null=None, title='', labels=None, ax=None, c=None, s=16, highlight_indices=None, highlight_c=None, highlight_label=None):
    """QQ-plot"""
    if labels is None:
        labels = ['', '']

    n = len(pval)
    x = -np.log10(np.arange(1,n+1)/(n+1))

    if ax is None:
        ax = setup_figure(4,4)

    ax.margins(x=0.02, y=0.05)
    args = {'s':s, 'edgecolor':'none', 'clip_on':False, 'alpha':1, 'rasterized':True}
    
    pval_idx_sorted = np.argsort(pval)
    log_pval_sorted = -np.log10(pval[pval_idx_sorted])

    ax.scatter(
        x,
        log_pval_sorted,
        c=c,
        zorder=30,
        label=labels[0],
        **args
    )
    
    if highlight_indices is not None:
        sorted_highlight_indices = [np.where(pval_idx_sorted == i)[0][0] for i in highlight_indices]
        ax.scatter(
            x[sorted_highlight_indices],
            log_pval_sorted[sorted_highlight_indices],
            c=highlight_c,
            zorder=40,
            label=highlight_label,
            **args
        )

    if pval_null is not None:
        assert len(pval)==len(pval_null)
        log_pval_sorted = -np.log10(np.sort(pval_null))
        ax.scatter(
            x,
            log_pval_sorted,
            c=[[0.5]*3],
            zorder=20,
            label=labels[1],
            **args
        )

    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True, min_n_ticks=5, nbins=4))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True, min_n_ticks=5, nbins=4))
    ax.set_xlabel('Expected -log$\mathregular{_{10}}$(p-value)', fontsize=14)
    ax.set_ylabel('Observed -log$\mathregular{_{10}}$(p-value)', fontsize=14)

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.set_xlim([0, xlim[1]])
    ax.set_ylim([0, ylim[1]])
    ci = 0.95
    xi = np.arange(1, n+1)
    clower = -np.log10(stats.beta.ppf((1-ci)/2, xi, xi[::-1]))
    cupper = -np.log10(stats.beta.ppf((1+ci)/2, xi, xi[::-1]))
    ax.fill_between(x, cupper, clower, color=[[0.8]*3], clip_on=True, rasterized=True)
    ax.plot([x[0], x[-1]], [x[0], x[-1]], '--', lw=1, color=[0.2]*3, zorder=50, clip_on=True, rasterized=True)
    #ax.spines['left'].set_position(('outward', 6))
    #ax.spines['bottom'].set_position(('outward', 6))
    ax.set_title('{}'.format(title), fontsize=12)
    if labels[0] != '':
        ax.legend(loc='upper left', fontsize=10, handlelength=0.5, handletextpad=0.33)
