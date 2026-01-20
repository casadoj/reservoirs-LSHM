import pandas as pd
from typing import Union, Dict, Optional
import matplotlib.pyplot as plt
from pathlib import Path


def plot_timeseries_MX(
    storage: pd.Series,
    elevation: pd.Series,
    outflow: Union[pd.Series, pd.DataFrame],
    meteo: Union[pd.Series, pd.DataFrame],
    max_storage: Dict,
    max_elevation: Dict,
    save: Optional[Union[str, Path]] = None,
    **kwargs
) -> None:
    """
    """
    
    lw = kwargs.get('lw', 1)
    figsize = kwargs.get('figsize', (12, 12))
    
    linestyles = [':', '--'] 
    fig, axes = plt.subplots(nrows=3, figsize=figsize, sharex=True)
    
    # storage
    ax1 = axes[0]
    ax1.plot(storage, c='C0', lw=lw, label='obs')
    for (key, value), ls in zip(max_storage.items(), linestyles):
        ax1.axhline(value, c='C0', ls=ls, lw=lw*.5, label=key)
    ax1.set_ylim(0, None)
    ax1.set_ylabel('storage (hm3)', color='C0')
    
    # elevation
    ax2 = ax1.twinx()
    ax2.plot(elevation, c='C1', lw=lw, label='obs')
    for (key, value), ls in zip(max_elevation.items(), linestyles):
        ax2.axhline(value, c='C1', ls=ls, lw=lw*.5, label=key)
    if 'zlim' in kwargs:
        ax2.set_ylim(kwargs['zlim'])
    ax2.set_ylabel('elevation (masl)', color='C1');
    
    # outflow
    ax = axes[1]
    outflow.plot(ax=ax, lw=lw)
    ax.set(ylabel='m3/s');
    
    # meteo
    ax = axes[2]
    meteo.plot(lw=1, ax=ax)
    ax.set(ylabel='mm');
    
    # title
    if 'title' in kwargs:
        axes[0].text(.5, 1.075, kwargs['title'],
                     ha='center', fontsize=12, transform=axes[0].transAxes);
    
    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight');
        plt.close()