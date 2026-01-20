import pandas as pd
from typing import Union, Dict, Optional
import matplotlib.pyplot as plt
from pathlib import Path


def plot_timeseries_AU(
    storage: pd.Series,
    elevation: pd.Series,
    max_storage: Dict,
    max_elevation: Dict,
    save: Optional[Union[str, Path]] = None,
    **kwargs
) -> None:
    """
    """
    
    lw = kwargs.get('lw', 1)
    figsize = kwargs.get('figsize', (16, 4))
    cmap = kwargs.get('cmap', 'plasma')
    c_V = kwargs.get('c1', '#1f77b4')
    c_Z = kwargs.get('c2', '#ff7f0e')
    
    linestyles = [':', '--'] 
    fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, figsize=figsize, gridspec_kw={'width_ratios': [1, 2.75], 'wspace': 0.1}, sharey=True)
    
    # Scatter plot (storage vs elevation)
    ax0.scatter(elevation, storage, c=storage.index, cmap=cmap, s=5, alpha=0.5)
    ax0.set_ylim(0, None)
    ax0.set_ylabel('storage (hm3)', color=c_V)
    ax0.set_xlabel('elevation (masl)', color=c_Z)
    for (key, value), ls in zip(max_storage.items(), linestyles):
        ax0.axhline(value, c=c_V, ls=ls, lw=lw*.5, label=key)
    for (key, value), ls in zip(max_elevation.items(), linestyles):
        ax0.axvline(value, c=c_Z, ls=ls, lw=lw*.5, label=key)
    ax0.tick_params(axis='x', colors=c_Z)
    ax0.tick_params(axis='y', colors=c_V)
    
    # storage
    ax1.plot(storage, c=c_V, lw=lw, label='obs')
    for (key, value), ls in zip(max_storage.items(), linestyles):
        ax1.axhline(value, c=c_V, ls=ls, lw=lw*.5, label=key)
    
    # elevation
    ax2 = ax1.twinx()
    ax2.plot(elevation, c=c_Z, lw=lw, label='obs')
    for (key, value), ls in zip(max_elevation.items(), linestyles):
        ax2.axhline(value, c=c_Z, ls=ls, lw=lw*.5, label=key)
    if 'zlim' in kwargs:
        ax2.set_ylim(kwargs['zlim'])
    ax2.set_ylabel('elevation (masl)', color=c_Z)
    ax2.tick_params(axis='y', colors=c_Z)
    
    # title
    if 'title' in kwargs:
        fig.text(.5, 0.925, kwargs['title'], ha='center', fontsize=12)
    
    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight');
        plt.close()