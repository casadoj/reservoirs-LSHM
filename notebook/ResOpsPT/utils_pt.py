import numpy as np
from numpy.polynomial.polynomial import Polynomial
import pandas as pd
from typing import Union, Dict, Optional
import matplotlib.pyplot as plt
from pathlib import Path


def fit_reservoir_curve(elevation: pd.Series, storage: pd.Series, deg=3, bias=.2):
    """
    """
    
    # merge time series, clean missing values and sort
    ts = pd.concat([elevation, storage], axis=1)
    ts.columns = ['elevation', 'storage']
    ts.dropna(how='any', inplace=True)
    ts.sort_values('elevation', inplace=True)

    # define training samples: remove duplicates in elevation
    curve = pd.Series(index=ts.elevation.unique(), name='storage', dtype=float)
    curve.index.name = 'elevation'
    for elev in curve.index:
        curve[elev] = ts.loc[ts.elevation == elev, ['storage']].median()[0]
    curve = curve.reset_index()
    
    # fit polynomial
    to_fit = True
    bias =.2
    while to_fit:

        # fit a degree-3 polynomial
        coefs = Polynomial.fit(curve.elevation, curve.storage, deg=3)

        # remove samples with a large error
        estimate = coefs(curve.elevation)
        error = estimate / curve.storage
        remove = (error < 1 - bias) | (error > 1 + bias)
        if remove.sum() > 0:
            # print(f'Remove {remove.sum()} points')
            curve = curve[~remove]
        else:
            to_fit = False
    
    
    # rmse = np.sqrt(np.mean((estimate - curve.storage)**2))
    # print('storage =', coefs)
    # print(f'RMSE = {rmse:.2f} hm3')
    
    return coefs


def plot_timeseries_PT(
    storage: Optional[pd.Series] = None,
    elevation: Optional[pd.Series] = None,
    outflow: Optional[pd.Series] = None,
    inflow: Optional[pd.Series] = None,
    max_storage: Optional[Dict] = None,
    max_elevation: Optional[Dict] = None,
    save: Optional[Union[str, Path]] = None,
    **kwargs
) -> None:
    """
    """
    
    lw = kwargs.get('lw', 1)
    s = kwargs.get('s', 5)
    figsize = kwargs.get('figsize', (16, 8))
    cmap = kwargs.get('cmap', 'plasma')
    c1 = kwargs.get('c1', '#1f77b4')
    c2 = kwargs.get('c2', '#ff7f0e')
    alpha = kwargs.get('alpha', .5)
    linestyles = [':', '--'] 

    df = pd.concat([serie for serie in [storage, elevation, inflow, outflow] if serie is not None], axis=1)
    start, end = df.first_valid_index(), df.last_valid_index()

    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=figsize, gridspec_kw={'width_ratios': [1, 3], 'wspace': 0.1})#, sharex=True, sharey=True)

    # scatter plot (storage vs elevation)
    if (storage is not None) and (elevation is not None):
        coefs = fit_reservoir_curve(elevation, storage, deg=3, bias=.2)
        elev_range = np.linspace(elevation.min(), elevation.max(), 100)
        ax[0,0].plot(elev_range, coefs(elev_range), lw=.5, ls=':', c='k', zorder=0)
        ax[0,0].scatter(elevation, storage, c='gray', marker='.', s=s, alpha=alpha)
        ax[0,0].set_ylabel('storage (hm3)', color=c1)
        ax[0,0].set_xlabel('elevation (masl)', color=c2)
        for (key, value), ls in zip(max_storage.items(), linestyles):
            ax[0,0].axhline(value, c=c1, ls=ls, lw=lw*.5, label=key)
        for (key, value), ls in zip(max_elevation.items(), linestyles):
            ax[0,0].axvline(value, c=c2, ls=ls, lw=lw*.5, label=key)
        ax[0,0].tick_params(axis='x', colors=c2)
        ax[0,0].tick_params(axis='y', colors=c1)
    else:
        ax[0,0].axis('off')

    # storage
    if storage is not None:
        ax[0,1].plot(storage, c=c1, lw=lw, label='obs')
        for (key, value), ls in zip(max_storage.items(), linestyles):
            ax[0,1].axhline(value, c=c1, ls=ls, lw=lw*.5, label=key)
        ax[0,1].set(
            xlim=(start, end),
            yticklabels=[]
        )
    else:
        ax[0,1].set_yticks([])

    # elevation
    ax2 = ax[0,1].twinx()
    ax2.plot(elevation, c=c2, lw=lw, label='obs')
    for (key, value), ls in zip(max_elevation.items(), linestyles):
        ax2.axhline(value, c=c2, ls=ls, lw=lw*.5, label=key)
    if 'zlim' in kwargs:
        ax2.set_ylim(kwargs['zlim'])
    ax2.set_ylabel('elevation (masl)', color=c2)
    ax2.tick_params(axis='y', colors=c2)
    ax2.set(
        xlim=(start, end)
    )

    # scatter plot (inflow vs outflow)
    if (inflow is not None) and (outflow is not None):
        ax[1,0].scatter(inflow, outflow, c='gray', marker='.', s=s, alpha=alpha)
        qmin = min(ax[1,0].get_xlim()[0], ax[1,0].get_ylim()[0])
        qmax = max(ax[1,0].get_xlim()[1], ax[1,0].get_ylim()[1])
        ax[1,0].plot([qmin, qmax], [qmin, qmax], lw=.5, ls=':', c='k', zorder=0)
        ax[1,0].set(xlim=(qmin, qmax), ylim=(qmin, qmax))
        ax[1,0].set_ylabel('outflow (m3/s)', color=c1)
        ax[1,0].set_xlabel('inflow (m3/s)', color=c2)
        ax[1,0].tick_params(axis='x', colors=c2)
        ax[1,0].tick_params(axis='y', colors=c1)
    else:
        ax[1,0].axis('off')

    # (in/out)flow
    if inflow is not None:
        ax[1,1].plot(inflow, c=c2, lw=lw, label='inflow')
    if outflow is not None:
        ax[1,1].plot(outflow, c=c1, lw=lw, label='outflow')
    ax[1,1].set(
        xlim=(start, end),
        yticklabels=[]
    )
    ax[1,1].legend(frameon=False)
    ax2 = ax[1,1].twinx()
    ax2.set(
        xlim=(start, end),
        ylim=ax[1,1].get_ylim(),
        ylabel='m3/s'
    )
    
    # title
    if 'title' in kwargs:
        fig.text(.5, 0.9, kwargs['title'], ha='center', fontsize=12)
    
    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight');
        plt.close()