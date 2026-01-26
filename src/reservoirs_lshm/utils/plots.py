import os
os.environ['USE_PYGEOS'] = '0'
import matplotlib as mpl
from matplotlib.axes import Axes
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.gridspec as gridspec
import geopandas as gpd
import seaborn as sns
import cartopy.crs as ccrs
import cartopy.feature as cf
import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path
from typing import Union, Dict, List, Tuple, Optional, Literal
from statsmodels.distributions.empirical_distribution import ECDF

# from utils import Decomposition
from .metrics import KGEmod, KGE, is_pareto_efficient


def plot_reservoir_map(
    geometry,
    volume: Optional[pd.Series] = None,
    area: Optional[pd.Series] = None,
    save: Union[str, Path] = None,
    **kwargs
):
    """Creates a maps where reservoirs are represented as dots. The size of the dots reflects the reservoir storage, and the colour the catchment area (if provided)
    
    Parameters:
    -----------
    geometry: gpd.GeoSeries
        Geometry of the points
    volume: pandas.Series
        Reservoir storage capacity (hm3)
    area: pandas.Series (optional)
        Reservoir catchment area (km2)
    save: str or pathlib.Path (optional)
        If provided, file where the plot will be saved
    """
    
    figsize = kwargs.get('figsize', (20, 5))
    title = kwargs.get('title', None)
    alpha = kwargs.get('alpha', .7)
    size = kwargs.get('size', 12)
    
    # set up plot
    fig, ax = plt.subplots(figsize=figsize, subplot_kw=dict(projection=ccrs.PlateCarree()))
    ax.add_feature(cf.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor='lightgray'), alpha=.5, zorder=0)
    ax.axis('off')
    
    # plot reservoir poitns
    s = np.sqrt(volume) if volume is not None else size
    c = np.sqrt(area) if area is not None else 'steelblue'
    cmap = 'coolwarm' if area is not None else None
    scatter = ax.scatter(
        geometry.x,
        geometry.y,
        s=s, #volume / 1000
        c=c, #, c=np.log10(icold.DOR_PC.replace(0, .1))
        cmap=cmap,
        alpha=alpha
    )
    # if area is not None:
    #     scatter = ax.scatter(
    #         geometry.x,
    #         geometry.y,
    #         s=np.sqrt(volume), #volume / 1000
    #         c=np.sqrt(area), #, c=np.log10(icold.DOR_PC.replace(0, .1))
    #         cmap='coolwarm',
    #         alpha=alpha
    #     )
    # else:
    #     scatter = ax.scatter(
    #         geometry.x,
    #         geometry.y,
    #         s=np.sqrt(volume), #volume / 1000
    #         color='steelblue',
    #         alpha=alpha
    #     )
    
    # title
    if title is not None:
        ax.text(.5, 1.125, title, horizontalalignment='center', verticalalignment='bottom', transform=ax.transAxes, fontsize=12)
    text = f'{len(geometry)} reservoirs'
    if volume is not None:
        text += '\n{0:.0f} km³'.format(volume.sum() / 1000)
    ax.text(.5, 1.02, text, horizontalalignment='center', verticalalignment='bottom', transform=ax.transAxes)
    
    # legend
    if area is not None:
        legend1 = ax.legend(*scatter.legend_elements(prop='colors', num=4, alpha=.5), title='catchment (km²)', bbox_to_anchor=[1.025, .6, .06, .25], frameon=False)
        ax.add_artist(legend1)
    if volume is not None:
        legend2 = ax.legend(*scatter.legend_elements(prop='sizes', num=4, alpha=.5), title='storage (km³)', bbox_to_anchor=[1.025, .1, .1, .5], frameon=False)
        ax.add_artist(legend2)

    # save
    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight')
        
        
def plot_attributes(
    df: pd.DataFrame,
    x: pd.Series,
    y: pd.Series,
    save: Optional[Union[Path, str]] = None,
    **kwargs
):
    """
    It creates maps (scatter plots) of the static attributes associated to specific points.

    Parameters:
    -----------
    df: pd.DataFrame
        table of attributes
    x: pd.Series
        coordinate X of the points in "df"
    y: pd.Series
        coordinate Y of the points in "df"
    save: optional, Path or str
        location where the plot will be saved. By default it is None and the plot won't be saved

    kwargs:
    -------
    figsize: List o Tuple
    ncols: int
    cmap: str
    alpha: float
    """

    # kwargs
    figsize = kwargs.get('figsize', (5, 4))
    ncols_max = kwargs.get('ncols', 3)
    cmap = kwargs.get('cmap', 'magma')
    alpha = kwargs.get('alpha', 1)
    s = kwargs.get('size', 5)
    extent = kwargs.get('extent', [-180, 180, -90, 90])
   
    proj = ccrs.PlateCarree()
    ncols, nrows = df.shape[1], 1
    if ncols > ncols_max:
        ncols, nrows = ncols_max, int(np.ceil(ncols / ncols_max))

    fig, axes = plt.subplots(ncols=ncols,
                             nrows=nrows,
                             figsize=(figsize[0] * ncols, figsize[1] * nrows),
                             subplot_kw={'projection': proj})
    for i, col in enumerate(df.columns):
        if nrows > 1:
            f, c = i // ncols, i % ncols
            ax = axes[f, c]
        else:
            c = i
            ax = axes[c]
        ax.add_feature(cf.NaturalEarthFeature('physical', 'land', '50m', edgecolor=None, facecolor='lightgray'), zorder=0)
        ax.set_extent(extent, crs=proj)
        sc = ax.scatter(x[df.index], y[df.index], cmap=cmap, c=df[col], s=s, alpha=alpha, label=col)
        cbar = plt.colorbar(sc, ax=ax, orientation='horizontal', shrink=.5)
        ax.set_title(' '.join(col.split('_')))
        ax.axis('off');
    
    if c < ncols - 1:
        for c_ in range(c + 1, ncols):
            axes[f, c_].axis('off')

    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight')
        
        
def compare_flows(
    storage: pd.Series,
    outflow: pd.Series,
    inflow1: pd.Series,
    inflow2: pd.Series,
    save: Union[str, Path] = None,
    **kwargs
):
    """
    It creates a figure with four plots. Three of them are scatter plots comparing pair-wisely the time series of outflow and inflows. The fourth plot shows the empirical cumulative of the three flow time series. Storage is only used in the colour scale.
    
    Parameters:
    -----------
    storage: pandas.Series
        Time series of observed reservoir storage
    outflow: pandas.Series
        Time series of observed reservoir outflow
    inflow1: pandas.Series
        One of the time series of reservoir inflow to be compared
    inflow2: pandas.Series
        The other time series of reservoir inflow to be compared
    save: string or pathlib.Path
        If provided, the file where the plot will be saved
        
    Returns:
    --------
    The plot is either printed on screen (default) or saved on disk (if "save" is provided)
    """
    
    figsize = kwargs.get('figsize', (8, 8))
    cmap = kwargs.get('cmap', 'coolwarm_r')
    s = kwargs.get('size', 4)
    a = kwargs.get('alpha', .5)
    scale = kwargs.get('scale', 'log')
    
    df = pd.concat((storage, outflow, inflow1, inflow2), axis=1)
    columns = df.columns
    
    fig, ax = plt.subplots(ncols=2, nrows=2, figsize=figsize, sharey=True)
    
    vmax = df.iloc[:,1:].max().max()
    r = len(str(int(vmax)))
    vmax = np.ceil(vmax / 10**r) * 10**r
    vmin = 0.001
    
    for c, x in enumerate([columns[2], columns[1]]):
        for r, y in enumerate([columns[1], columns[3]]):
            if c == 1 and r == 0:
                continue
                
            sct = ax[r,c].scatter(df[x], df[y], c=df[columns[0]], cmap=cmap, s=s, alpha=a)
            ax[r,c].plot([vmin, vmax], [vmin, vmax], '--k', lw=.5, zorder=0)
            
            try:
                kge, alpha, beta, rho = KGE(df[x], df[y])
                ax[r,c].text(.01, .99, f'KGE={kge:.2f}', transform=ax[r,c].transAxes, ha='left', va='top')
                ax[r,c].text(.01, .94, f'α={alpha:.2f}', transform=ax[r,c].transAxes, ha='left', va='top')
                ax[r,c].text(.01, .89, f'β={beta:.2f}', transform=ax[r,c].transAxes, ha='left', va='top')
                ax[r,c].text(.01, .84, f'ρ={rho:.2f}', transform=ax[r,c].transAxes, ha='left', va='top')
            except:
                pass
            
            ax[r,c].set(xlim=(vmin, vmax),
                        xscale=scale, 
                        ylim=(vmin, vmax),
                        yscale=scale)
            if r == 1:
                ax[r,c].set_xlabel('{0} (m3/s)'.format(x))
            if c == 0:
                ax[r,c].set_ylabel('{0} (m3/s)'.format(y))
                if r == 0:
                    ax[r,c].set_xticklabels([])
    
    # ecdf
    for col, ls in zip(columns[1:], ['-', '--', ':']):
        ecdf = ECDF(df[col])
        ax[0,1].plot(ecdf.y, ecdf.x, c='k', ls=ls, lw=1, label=col)
    ax[0,1].set(xlabel='ECDF (-)',
                ylabel='flow (m3/s)',
                xlim=(-.02, 1.02),
                ylim=(vmin, vmax),
                yscale=scale)
    ax[0,1].legend(frameon=False);

    cbar_ax = fig.add_axes([0.33, 0.02, 0.33, 0.015])
    plt.colorbar(sct, cax=cbar_ax, orientation='horizontal', label='storage (hm3)')
    
    if 'title' in kwargs:
        fig.text(0.5, .925, kwargs['title'], ha='center', va='top', fontsize=12)
    
    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight')
        plt.close(fig)
        
        
def create_cmap(
    cmap: str,
    bounds: List,
    name: str = '',
    specify_color: Tuple = None
):
    """Given the name of a colour map and the boundaries, it creates a discrete colour ramp for future plots
    
    Inputs:
    ------
    cmap:          string. Matplotlib's name of a colourmap. E.g. 'coolwarm', 'Blues'...
    bounds:        list. Values that define the limits of the discrete colour ramp
    name:          string. Optional. Name given to the colour ramp
    specify_color: tuple (position, color). It defines a specific color for a specific position in the colour scale. Position must be an integer, and color must be either a colour name or a tuple of 4 floats (red, gren, blue, transparency)
    
    Outputs:
    --------
    cmap:   List of colours
    norm:   List of boundaries
    """
    
    cmap = plt.get_cmap(cmap)
    cmaplist = [cmap(i) for i in range(cmap.N)]
    if specify_color is not None:
        cmaplist[specify_color[0]] = specify_color[1]
    cmap = mpl.colors.LinearSegmentedColormap.from_list(name, cmaplist, cmap.N)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    
    return cmap, norm

def plot_resops(
    storage: pd.Series = None,
    elevation: pd.Series = None,
    inflow: pd.Series = None,
    outflow: pd.Series = None,
    capacity: Union[List[float], float] = None,
    level: Union[List[float], float] = None,
    save: Union[str, Path] = None,
    **kwargs
):
    """It creates a plot with two graphs that shows the reservoir time series. The first graph is the storage-elevation curve of the reservoir (if both storage and elevation time series are available). The second graph is the time series of storage, inflow and outflow.
    
    Parameters:
    -----------
    storage: pd.Series
        Time series of reservoir storage (hm3)
    elevation: pd.Series
        Time series of reservoir level (masl)
    inflow: pd.Series
        Time series of reservoir inflow (m3/s)
    outflow: pd.Series
        Time series of reservoir outflow (m3/s)
    capacity: Union[List[float], float]
        Values of storage capacity that will be plotted as horizontal lines in the plots
    level: float or list of floats
        Values of reservoir level that will be plotted as vertical lies in the scatter plot
    save: Union[str, Path]
        File name where the plot will be saved. By default is None and the plot is not saved.
    kwargs:
        figsize: tuple
            Size of the plot
        xlim: tuple
            Limits of the X axis (time) in the time series plot
        ylim: tuple
            Limites of the Y axis (storage) in both plots. They share the Y axis
        title: str
            If given, title of the plot
    """
                
    # Create the figure and define the grid specification
    fig = plt.figure(figsize=kwargs.get('figsize', (20, 6)))
    gs = gridspec.GridSpec(1, 3)

    # Create the first graph in the left-most third
    ax1 = plt.subplot(gs[0])
    if isinstance(elevation, pd.Series) and isinstance(storage, pd.Series):
        ax1.scatter(elevation, storage, s=1, c=elevation.index, cmap='Greys')
    if capacity is not None:
        if isinstance(capacity, float):
            capacity = [capacity]
        for ls, value in zip(['-', '--', ':'], capacity):
            ax1.axhline(value, lw=.8, c='k', ls=ls)
    if level is not None:
        if isinstance(elevation, float):
            level = [level]
        for ls, value in zip(['-', '--', ':'], level):
            ax1.axvline(value, lw=.8, c='k', ls=ls)
    ax1.set(xlabel='elevation (m)',
            ylabel='storage (hm3)')

    # Create the second graph in the remaining area
    ax2 = plt.subplot(gs[1:], sharey=ax1)
    if isinstance(storage, pd.Series):
        ax2.fill_between(storage.index, storage, color='lightgray', alpha=.5, label='storage')
    if capacity is not None:
        for ls, value in zip(['-', '--', ':'], capacity):
            ax2.axhline(value, lw=.8, c='k', ls=ls)
    ax2.set(xlim=kwargs.get('xlim', (storage.index[0], storage.index[-1])),
            ylim=(0, None))

    ax2_ = ax2.twinx()
    if isinstance(inflow, pd.Series):
        ax2_.plot(inflow, lw='.8', c='steelblue', label='inflow')
    if isinstance(outflow, pd.Series):
        ax2_.plot(outflow, lw='.8', c='indianred', label='outflow')
    ax2_.set(ylim=(0, None),
             ylabel='flow (m3/s)');
    
    fig.legend(ncol=3, loc=8, frameon=False, bbox_to_anchor=[.5, .0, .3, .1])
    fig.text(.5, .925, kwargs.get('title', None), fontsize=15, horizontalalignment='center')
    
    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight')
        plt.close(fig)
        
        
# def plot_decomposition(obs: Decomposition, sim: Decomposition, id: int, lims: List[float] = [.1, .67, .97], save: Union[str, Path] = None, **kwargs):
#     """It creates a figure that compares the decomposition of the observed and simulated time series. The figure is composed of 4 plots: original time series, trend, seasonality and residuals. Each plot includes the performance in terms of modified KGE and its components.
    
#     Parameters:
#     -----------
#     obs:       Decomposition
#         The result of decomposing the observed timeseries with the function utils.decompose_timeseries
#     sim:       Decomposition
#         The result of decomposing the simulated timeseries with the function utils.decompose_timeseries
#     id:        int
#         Identification number of the station to be plotted
#     lims:      List[float]
#         Values of the conservative (clim), normal (nlim) and flood (flim) limits in the LISFLOOD parameterization. If not provided, it takes the default values in GloFAS
#     save:      bool
#         Whether to save of not the figure. By default is None and the figure is note saved
        
#     kwargs:
#     -------
#     figsize:   Tuple (2,)
#         Size of each of the individual plots
#     lw:        float
#         Line width
#     title:     str
#         If provided, the title of the figure
#     """
    
#     # kwargs
#     figsize = kwargs.get('figsize', (3, 6))
#     lw = kwargs.get('lw', 1)
    
#     # setup
#     bbox_props = dict(boxstyle='round, pad=0.05', facecolor='w', edgecolor='none', alpha=.666)
#     ncols = 4
#     fig, axes = plt.subplots(figsize=(figsize[0] * ncols, figsize[1]), ncols=ncols, tight_layout=True, sharey=True)#, sharex=True)
    
#     components = [attr for attr in dir(obs) if not attr.startswith('_')]
#     for ax, comp in zip(axes, ['original', 'trend', 'seasonal', 'residual']):
#         # extract series
#         obs_comp = getattr(obs, comp)[id]
#         sim_comp = getattr(sim, comp)[id]
        
#         # line plots
#         ax.plot(obs_comp, obs_comp.index, lw=lw, label='obs')
#         ax.plot(sim_comp, sim_comp.index, lw=lw, label='sim')
        
#         # add performance as text
#         performance = KGEmod(obs_comp, sim_comp)
#         performance = ['-∞' if x < -10 else '∞' if x > 10 else np.round(x, 2) for x in performance]
#         ax.text(.02, .99, "KGE' = {0}\nα = {1}\nβ = {2}\nρ = {3}".format(*performance),
#                 va='top', transform=ax.transAxes, bbox=bbox_props)
        
#         # settings
#         if comp in ['original', 'trend']:
#             for lim in lims:
#                 ax.axvline(lim, ls=':', c='k', lw=.5)
#             ax.set(xlim=(-.02, 1.02),
#                    ylim=(obs_comp.index.max(), obs_comp.index.min()))
#         else:
#             ax.set(xlim=(-.51, .51))
#         ax.set(xlabel='Snorm',
#                title=comp);
        
#     fig.legend(*ax.get_legend_handles_labels(), ncol=2, loc=8, frameon=False, bbox_to_anchor=[.4, -.025, .2, .1])
#     if 'title' in kwargs:
#         fig.text(.5, 1, kwargs['title'], fontsize=11, ha='center');

#     if save is not None:
#         plt.savefig(save, dpi=300, bbox_inches='tight');


def plot_decomposition(
    sim: pd.DataFrame,
    obs: pd.DataFrame = None,
    lims: List[float] = None,
    save: Union[str, Path] = None,
    **kwargs
):
    """It creates a figure that compares the decomposition of the observed and simulated time series. The figure is composed of 4 plots: original time series, trend, seasonality and residuals. Each plot includes the performance in terms of modified KGE and its components.
    
    Parameters:
    -----------
    sim:       Decomposition
        The result of decomposing the simulated timeseries with the function utils.decompose_timeseries
    obs:       Decomposition
        The result of decomposing the observed timeseries with the function utils.decompose_timeseries
    lims:      List[float]
        Values of the conservative (clim), normal (nlim) and flood (flim) limits in the LISFLOOD parameterization. If not provided, it takes the default values in GloFAS
    save:      bool
        Whether to save of not the figure. By default is None and the figure is note saved
        
    kwargs:
    -------
    figsize:   Tuple (2,)
        Size of each of the individual plots
    lw:        float
        Line width
    title:     str
        If provided, the title of the figure
    """
    
    # kwargs
    figsize = kwargs.get('figsize', (3, 6))
    lw = kwargs.get('lw', 1)
    xlim = kwargs.get('xlim', [(-.02, 1.02), (-.51, .51)])
    xlabel = kwargs.get('xlabel', 'Snorm')
    c = kwargs.get('color', ['C1', 'C0'])
    
    # setup
    bbox_props = dict(boxstyle='round, pad=0.05', facecolor='w', edgecolor='none', alpha=.666)
    ncols = 4
    fig, axes = plt.subplots(figsize=(figsize[0] * ncols, figsize[1]), ncols=ncols, tight_layout=True, sharey=True)#, sharex=True)
    
    components = [attr for attr in dir(obs) if not attr.startswith('_')]
    for ax, comp in zip(axes, sim.columns):
        if obs is not None:
            # observed time series
            obs_comp = obs[comp] #getattr(obs, comp)[id]
            ax.plot(obs_comp, obs_comp.index, c=c[1], lw=lw, label='obs')
            
        # simulated time series
        sim_comp = sim[comp] #getattr(sim, comp)[id]
        ax.plot(sim_comp, sim_comp.index, c=c[0], lw=lw, label='sim')
            
        if obs is not None:        
            # add performance as text
            performance = KGEmod(obs_comp, sim_comp)
            performance = ['-∞' if x < -10 else '∞' if x > 10 else np.round(x, 2) for x in performance]
            ax.text(.02, .99, "KGE' = {0}\nα = {1}\nβ = {2}\nρ = {3}".format(*performance),
                    va='top', transform=ax.transAxes, bbox=bbox_props)
        
        # settings
        if comp in ['original', 'trend']:
            if lims is not None:
                for lim in lims:
                    ax.axvline(lim, ls=':', c='k', lw=.5)
            ax.set(xlim=xlim[0],
                   ylim=(sim_comp.index.max(), sim_comp.index.min()))
        else:
            ax.set(xlim=xlim[1])
        ax.set(xlabel=xlabel,
               title=comp);
        
    fig.legend(*ax.get_legend_handles_labels(), ncol=2, loc=8, frameon=False, bbox_to_anchor=[.4, -.025, .2, .1])
    if 'title' in kwargs:
        fig.text(.5, 1, kwargs['title'], fontsize=11, ha='center');

    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight')
        plt.close(fig)
        

# def storage_outflow(storage: pd.Series, outflow: pd.Series, storage2: pd.Series = None, outflow2: pd.Series = None, s_lims: List = None, q_lims: List = None, save: Union[Path, str] = None, **kwargs):
#     """It creates a figure that compares the storage and outflow time series. The figure is composed of three plots. In the center, a scatter plot of storage versus outflow; if the storage and outflow limits are provided, a line
#     represents the reference LISFLOOD routine. On top, a plot shows the density function (kernel density estimation) of storage. On the right, a plot shows the density function (kernel density estimation) of outflow.
    
#     Parameters:
#     -----------
#     storage:   pd.Series
#         Series of reservoir storage. By default, it should be storage relative to the total reservoir capacity. It is supposed to be simulated storage
#     outflow:   pd.Series
#         Series of reservoir outflow. By default, it should be relative to the non-damaging outflow. It is supposed to be simulated outflow
#     storage2:  pd.Series
#         Series of reservoir storage. By default, it should be storage relative to the total reservoir capacity. It is supposed to be observed storage
#     outflow2:  pd.Series
#         Series of reservoir outflow. By default, it should be relative to the non-damaging outflow. It is supposed to be observed outflow
#     s_lims:    List
#         Storage limits (conservative, 2 times conservation, normal, adjusted normal, flood) used in the LISFLOOD reservoir routine
#     q_lims:    List
#         Outflow limits (minimum, minimum, normal adjusted, normal adjusted, non-damaging) used in the LISFLOOD reservoir routine
#     save:      Union[str, Path]
#         Path where to save the figure
    
#     Keyword arguments:
#     ------------------
#     alpha:     float
#         Transparency in the scatter plot
#     color:     list(2,)
#         Colours to be used in the simulated and observed data
#     size:      float
#         Point size in the scatter plot
#     title:     str
#         Title of the figure
#     xlabel:    str
#         Label of the X axis in the scatter plot
#     ylabel:    str
#         Label of the Y axis in the scatter plot
#     ymin:      float
#         Minimum value of the Y axis in the scatter plot
#     """
    
#     # extract kwargs
#     s = kwargs.get('size', .5)
#     a = kwargs.get('alpha', .05)
#     c = kwargs.get('color', ['C1', 'C0'])
#     ymin = kwargs.get('ymin', -.1)
#     xlabel = kwargs.get('xlabel', 'relative storage (-)')
#     ylabel = kwargs.get('ylabel', 'relative outflow (-)')
    
#     if storage2 is not None:
#         if storage2.isnull().all():
#             storage2 = None
#     if outflow2 is not None:
#         if outflow2.isnull().all():
#             outflow2 = None
    
#     if (s_lims is not None) & (q_lims is not None):
#         assert len(s_lims) == len(q_lims), 'The length of "s_lims" and "q_lims" must be the same.'
    
#     # Create the figure and set the size
#     fig = plt.figure(figsize=(5, 5))
#     gs = gridspec.GridSpec(2, 2, height_ratios=[1, 4], width_ratios=[4, 1])
#     if 'title' in kwargs:
#         fig.text(.95, .95, kwargs['title'], ha='right', va='top')  

#     # scatter plot outflow vs storage
#     ax1 = plt.subplot(gs[1, 0])
#     ax1.scatter(storage, outflow, c=c[0], s=s, alpha=a)
#     ax1.scatter(-1, -1, c=c[0], s=s, label=kwargs.get('label1', 'sim'))
#     if (storage2 is not None) & (outflow2 is not None):
#         ax1.scatter(storage2, outflow2, c=c[1], s=s, alpha=a)
#         ax1.scatter(-1, -1, s=s, c=c[1], label=kwargs.get('label2', 'obs'))
#     if (s_lims is not None) & (q_lims is not None):
#         ax1.plot(s_lims, q_lims, c='k', lw=1, zorder=0, label='routine');
#         for s, q in zip(s_lims, q_lims):
#             ax1.hlines(q, xmin=-.02, xmax=s, color='k', ls=':', lw=.5, zorder=10)
#             ax1.vlines(s, ymin=ymin, ymax=q, color='k', ls=':', lw=.5, zorder=10)
#     ax1.set(ylim= (ymin, None),
#             xlim=(-.02, 1.02), 
#             xlabel=xlabel,
#             ylabel=ylabel)
#     ax1.spines[['top', 'right']].set_visible(False)
#     ax1.legend(frameon=False, loc=2)
    
#     # densidy distribution: storage
#     ax2 = plt.subplot(gs[0, 0])
#     sns.kdeplot(storage, color=c[0], fill=True, ax=ax2)
#     if storage2 is not None:
#         sns.kdeplot(storage2, color=c[1], fill=True, ax=ax2)
#         kge = KGEmod(storage, storage2)[0]
#         ax2.text(.5, -.2, f"KGE' = {kge:.2f}", fontsize=9, va='center', ha='center', transform=ax2.transAxes)
#     if s_lims is not None:
#         for s in s_lims:
#             ax2.axvline(s, color='k', ls=':', lw=.5, zorder=10)
#     ax2.spines[['top', 'left', 'right']].set_visible(False)
#     ax2.set(xlim=(-.02, 1.02),
#             ylabel=None,
#             xlabel=None)
#     ax2.set_yticks([])
#     ax2.set_xticklabels([])

#     # density distribution: outflow
#     ax3 = plt.subplot(gs[1, 1])
#     sns.kdeplot(y=outflow, color=c[0], fill=True, ax=ax3)
#     if outflow2 is not None:
#         sns.kdeplot(y=outflow2, color=c[1], fill=True, ax=ax3)
#         kge = KGEmod(outflow, outflow2)[0]
#         ax3.text(-0.2, .55, f"KGE' = {kge:.2f}", fontsize=9, va='center', ha='center', rotation=90, transform=ax3.transAxes)
#     if q_lims is not None:
#         for q in q_lims:
#             ax3.axhline(q, color='k', ls=':', lw=.5, zorder=10)
#     ax3.spines[['top', 'bottom', 'right']].set_visible(False)
#     ax3.set(ylim=(ymin, None),
#             xlabel=None,
#             ylabel=None)
#     ax3.set_xticks([])
#     ax3.set_yticklabels([])

#     # Adjust the spacing between subplots
#     fig.tight_layout();
    
#     if save is not None:
#         plt.savefig(save, dpi=300, bbox_inches='tight')
        
        
def reservoir_scatter(
    sim: pd.DataFrame, 
    x: str, 
    y: str, 
    obs: Optional[pd.DataFrame] = None, 
    x_thr: Optional[List] = None, 
    y_thr: Optional[List] = None, 
    legend: bool = True, 
    ax: Optional[Axes] = None, 
    **kwargs
):
    """It creates a figure that compares the storage and outflow time series. The figure is composed of three plots. In the center, a scatter plot of storage versus outflow; if the storage and outflow limits are provided, a line represents the reference LISFLOOD routine. On top, a plot shows the density function (kernel density estimation) of storage. On the right, a plot shows the density function (kernel density estimation) of outflow.
    
    Parameters:
    -----------
    sim:   pd.DataFrame
        Simulated time series of reservoir behaviour. It should contain colums "x" and "y"
    x:     str
        Column of "sim" (and "obs") to be used in the X axis
    y:     str
        Column of "sim" (and "obs") to be used in the Y axis
    obs:   pd.DataFrame (optional)
        Oberved time series of reservoir behaviour. It should contain colums "x" and "y"
    x_thr:    List (optional)
        Thresholds in the LISFLOOD reservoir routine to be used in the X axis
    y_thr:    List (optional)
        Thresholds in the LISFLOOD reservoir routine to be used in the Y axis
    legend:   bool
        Whether to plot the legend or not
    ax:       Axes (optional)
        Matplotlib axes in which to insert the plot
    
    Keyword arguments:
    ------------------
    alpha:     float
        Transparency in the scatter plot
    color:     List(2,)
        Colours to be used in the simulated and observed data
    figsize    Tuple(2,)
        Size of the figure
    labels:    List(2,)
        Labels of the datasets "sim" and "obs"
    size:      float
        Point size in the scatter plot
    xlabel:    str
        Label of the X axis in the scatter plot
    xlim:      Tuple(2,)
        Limits of the X axis
    xticklabels: bool
        Whether to include values in the X ticks or not
    ylabel:    str
        Label of the Y axis in the scatter plot
    ylim:      Tuple(2,)
        Limits of the Y axis
    yticklabels: bool
        Whether to include values in the Y ticks or not
    """
    
    # extract kwargs
    a = kwargs.get('alpha', .05)
    c = kwargs.get('color', ['C1', 'C0'])
    figsize = kwargs.get('figsize', (4, 4))
    labels = kwargs.get('labels', ['sim', 'obs'])
    s = kwargs.get('size', .5)
    xlabel = kwargs.get('xlabel', x)
    xlim = kwargs.get('xlim', (-.1, None))
    xticklabels = kwargs.get('xticklabels', True)
    ylabel = kwargs.get('ylabel', y)
    ylim = kwargs.get('ylim', (-.1, None))
    yticklabels = kwargs.get('yticklabels', True)
    
    if (x_thr is not None) & (y_thr is not None):
        assert len(x_thr) == len(y_thr), 'The length of "x_thr" and "y_thr" must be the same.'
    
    # Create the figure and set the size
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # scatter plot outflow vs storage
    ax.scatter(sim[x], sim[y], c=c[0], s=s, alpha=a, zorder=2)
    ax.scatter(-1, -1, c=c[0], s=s, label=labels[0])
    if obs is not None:
        ax.scatter(obs[x], obs[y], c=c[1], s=s, alpha=a, zorder=1)
        ax.scatter(-1, -1, s=s, c=c[1], label=labels[1])
    if (x_thr is not None) & (y_thr is not None):
        #ax.plot(x_thr, y_thr, c='k', lw=1, zorder=0, label='routine');
        for s, q in zip(x_thr, y_thr):
            ax.hlines(q, xmin=-.02, xmax=s, color='k', ls=':', lw=.5, zorder=10)
            ax.vlines(s, ymin=ylim[0], ymax=q, color='k', ls=':', lw=.5, zorder=10)
    ax.set(xlim=xlim,
           ylim=ylim,
           xlabel=xlabel,
           ylabel=ylabel)
    if xticklabels is False:
        ax.set_xticklabels([])
    if yticklabels is False:
        ax.set_yticklabels([])
    ax.spines[['top', 'right']].set_visible(False)
    
    if legend:
        ax.legend(frameon=False, loc=2)
 
        
def reservoir_kde(
    sim: pd.DataFrame, 
    obs: pd.DataFrame = None, 
    x: str = None, 
    y: str = None, 
    thr: List = None, 
    ax: Axes = None, 
    **kwargs
):
    """It creates a figure that compares the storage and outflow time series. The figure is composed of three plots. In the center, a scatter plot of storage versus outflow; if the storage and outflow limits are provided, a line represents the reference LISFLOOD routine. On top, a plot shows the density function (kernel density estimation) of storage. On the right, a plot shows the density function (kernel density estimation) of outflow.
    
    Parameters:
    -----------
    sim:   pd.DataFrame
        Simulated time series of reservoir behaviour. It should contain colums "x" and "y"
    obs:   pd.DataFrame
        Oberved time series of reservoir behaviour. It should contain colums "x" and "y"
    x:     str
        Column of "sim" (and "obs") for which the distribution will be computed. The plot will be horizontal
    y:     str
        Column of "sim" (and "obs") for which the distribution will be computed. The plot will be vertical
    thr:    List
        Thresholds in the LISFLOOD reservoir routine
    ax:       Axes
        Matplotlib axes in which to insert the plot
    
    Keyword arguments:
    ------------------
    color:     List(2,)
        Colours to be used in the simulated and observed data
    figsize    Tuple(2,)
        Size of the figure
    xlabel:    str
        Label of the X axis in the scatter plot
    xlim:      Tuple(2,)
        Limits of the X axis
    xticklabels: bool
        Whether to include values in the X ticks or not
    ylabel:    str
        Label of the Y axis in the scatter plot
    ylim:      Tuple(2,)
        Limits of the Y axis
    yticklabels: bool
        Whether to include values in the Y ticks or not
    """
    
    assert (x != y) & ((x is not None) | (y is not None)), 'Either "x" or "y" must indicate a column in "sim".'
    
    # extract kwargs
    c = kwargs.get('color', ['C1', 'C0'])
    figsize = kwargs.get('figsize', (4, 4))
    xlim = kwargs.get('xlim', (-.1, None))
    xlabel = kwargs.get('xlabel', None)
    xticklabels = kwargs.get('xticklabels', True)
    ylabel = kwargs.get('ylabel', None)
    ylim = kwargs.get('ylim', (-.1, None))
    yticklabels = kwargs.get('yticklabels', True)
    
    # Create the figure and set the size
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)  
    
    # density distribution
    if x is not None:
        sns.kdeplot(sim[x], color=c[0], fill=True, ax=ax)
        if obs is not None:
            if not obs[x].isnull().all():
                sns.kdeplot(obs[x], color=c[1], fill=True, ax=ax)
                kge = KGEmod(obs[x], sim[x])[0]
                ax.text(.5, -.2, f"KGE' = {kge:.2f}", fontsize=9, va='center', ha='center', transform=ax.transAxes)
        if thr is not None:
            for x in thr:
                ax.axvline(x, color='k', ls=':', lw=.5, zorder=10)
        ax.spines[['top', 'left', 'right']].set_visible(False)
        ax.set(xlim=xlim,
                 ylabel=None,
                 xlabel=None)
        ax.set_yticks([])
        if xticklabels is False:
            ax.set_xticklabels([])
    elif y is not None:
        sns.kdeplot(y=sim[y], color=c[0], fill=True, ax=ax)
        if obs is not None:
            if not obs[y].isnull().all():
                sns.kdeplot(y=obs[y], color=c[1], fill=True, ax=ax)
                kge = KGEmod(obs[y], sim[y])[0]           
                ax.text(-0.2, .55, f"KGE' = {kge:.2f}", fontsize=9, va='center', ha='center', rotation=90, transform=ax.transAxes)
        if thr is not None:
            for y in thr:
                ax.axhline(y, color='k', ls=':', lw=.5, zorder=10)
        ax.spines[['top', 'bottom', 'right']].set_visible(False)
        ax.set(ylim=ylim,
               xlabel=None,
               ylabel=None)
        ax.set_xticks([])
        if yticklabels is False:
            ax.set_yticklabels([])

            
def reservoir_analysis(
    sim: pd.DataFrame,
    obs: Optional[pd.DataFrame] = None,
    x1: str = 'storage',
    x2: str = 'inflow',
    y: str = 'outflow',
    x_thr: Optional[List] = None,
    y_thr: Optional[List] = None,
    save: Optional[Union[Path, str]] = None,
    **kwargs
):
    """It creates a figure that compares the storage and outflow time series. The figure is composed of three plots. In the center, a scatter plot of storage versus outflow; if the storage and outflow limits are provided, a line
    represents the reference LISFLOOD routine. On top, a plot shows the density function (kernel density estimation) of storage. On the right, a plot shows the density function (kernel density estimation) of outflow.
    
    Parameters:
    -----------
    sim:   pd.DataFrame
        Simulated time series of reservoir behaviour. It should contain colums "x" and "y"
    obs:   pd.DataFrame (optional)
        Oberved time series of reservoir behaviour. It should contain colums "x" and "y"
    x1:     str
        Column of "sim" (and "obs") that will be used in the X axis of the first scatter plot
    x2:     str
        Column of "sim" (and "obs") that will be used in the X axis of the second scatter plot
    y:     str
        Column of "sim" (and "obs") that will be used in the Y axis of both scatter plots
    x_thr:    List (optional)
        Thresholds in the LISFLOOD reservoir routine to be used in the "x1" variable
    y_thr:    List (optional)
        Thresholds in the LISFLOOD reservoir routine to be used in the "y" axis
    save:      string or pathlib.Path (optional)
        Path where to save the figure
    
    Keyword arguments:
    ------------------
    alpha:     float
        Transparency in the scatter plot
    color:     list(2,)
        Colours to be used in the simulated and observed data
    figsize:   Tuple(2,)
        Size of the figure
    labels:    List(2,)
        Labels of the datasets "sim" and "obs"
    size:      float
        Point size in the scatter plot
    title:     str
        Title of the figure
    x1lim:     Tuple(2,)
        Limits of the "x1" axis
    """
    
    # extract kwargs
    a = kwargs.get('alpha', .05)
    c = kwargs.get('color', ['C1', 'C0'])
    figsize = kwargs.get('figsize', (9, 5))
    labels = kwargs.get('labels', ['sim', 'obs'])
    s = kwargs.get('size', .5)
    x1_lim = kwargs.get('x1lim', (-.02, 1.02))
    
    if (x_thr is not None) & (y_thr is not None):
        assert len(x_thr) == len(y_thr), 'The length of "x_thr" and "y_thr" must be the same.'
    
    # Create the figure and set the size
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(2, 3, height_ratios=[1, 4], width_ratios=[4, 4, 1])
    if 'title' in kwargs:
        fig.text(.95, .95, kwargs['title'], ha='right', va='top')  

    # scatter plot: x1 vs y
    ax10 = plt.subplot(gs[1, 0])
    reservoir_scatter(
        sim, 
        x1, 
        y, 
        obs, 
        x_thr=x_thr, 
        y_thr=y_thr, 
        xlim=x1_lim, 
        ax=ax10, 
        legend=False,
        size=s, 
        alpha=a, 
        color=c, 
        labels=labels
    )
    
    # scatter plot: x2 vs y
    ax11 = plt.subplot(gs[1, 1])
    reservoir_scatter(
        sim, 
        x2,
        y, 
        obs, 
        x_thr=y_thr, 
        y_thr=y_thr, 
        ax=ax11, 
        legend=False, 
        ylim=ax10.get_ylim(), 
        ylabel='', 
        yticklabels=False,
        size=s, 
        alpha=a, 
        color=c, 
        labels=labels
    )
    ax11.plot(ax10.get_ylim(), ax10.get_ylim(), c='k', lw=.5, ls=':', zorder=0)
    
    # density distribution: x1
    ax00 = plt.subplot(gs[0, 0])
    reservoir_kde(
        sim, 
        obs, 
        x=x1, 
        thr=x_thr, 
        ax=ax00, 
        xlim=ax10.get_xlim(), 
        xticklabels=False,   
        color=c
    )
    
    # density distribution: x2
    ax01 = plt.subplot(gs[0, 1])
    reservoir_kde(
        sim, 
        obs, 
        x=x2, 
        thr=y_thr, 
        ax=ax01, 
        xlim=ax11.get_xlim(), 
        xticklabels=False,
        color=c
    )

    # density distribution: y
    ax12 = plt.subplot(gs[1, 2])
    reservoir_kde(
        sim, 
        obs, 
        y=y, 
        thr=y_thr, 
        ax=ax12, 
        ylim=ax10.get_ylim(), 
        yticklabels=False, 
        color=c
    )
    
    fig.legend(*ax10.get_legend_handles_labels(), frameon=False, ncol=3, loc=8, bbox_to_anchor=[.25, -.04, .5, .05])
    
    # Adjust the spacing between subplots
    fig.tight_layout();
    
    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight')
        plt.close(fig)
        
        
def maps_performance(
    x: pd.Series, 
    y: pd.Series, 
    performance: pd.DataFrame, 
    s: Union[pd.Series, int] = None, 
    polygons: gpd.GeoDataFrame = None, 
    save: Union[Path, str] = None, 
    **kwargs
):
    """It creates a figure that contains 4 maps with the KGE and its 3 components.
    
    Inputs:
    -------
    x: pd.Series
        X coordinate of the points 
    y: pd.Series
        Y coordinate of the points
    performance: pd.DataFrame
        Performance of each of the points. It must contain at least 4 columns: 'KGE', 'r', 'alpha', 'beta'
    s: Union[pd.Series, int]
        Size of each point. In case of an integer all points will have the same size. In case of a pd.Series every point will have an specific size
    save: Union[Path, str]
        File in which the plot will be saved
        
    kwargs:
    -------
    figsize: Tuple
        Size of the figure
    proj: 
        Projection of the map
    extent: List
        Extension of the map [xmin, xmax, ymin, ymax]
    alpha: float
        Transparency of the points
    title: str
        Title of the figure
    """
    
    figsize = kwargs.get('figsize', (15, 7))
    proj = kwargs.get('proj', ccrs.PlateCarree())
    extent = kwargs.get('extent', None)#[-125, -67.5, 24, 51]
    alpha = kwargs.get('alpha', .8)
    if s is None:
        s = 5
    
    fig, axes = plt.subplots(ncols=2, nrows=2, figsize=figsize, tight_layout=True, subplot_kw={'projection': proj})

    for i, metric in enumerate(['KGE', 'r', 'alpha', 'beta']):

        if metric == 'KGE':
            cmap, norm = create_cmap('Spectral', [-100, -1, -.75, -.5, -.25 ,0, .25, .5, .75, 1])
        elif metric in ['alpha', 'beta']:
            cmap, norm = create_cmap('RdBu', [1e-6, 1/16, 1/8, 1/4, 1/2, 1/1.2, 1.2, 2, 4, 8, 16, 1e6])
        elif metric == 'r':
            cmap, norm = create_cmap('Spectral', [-1, -.75, -.5, -.25, 0, .25, .5, .75, 1])

        # background map
        ax = axes[int(i / 2), i % 2]
        ax.add_feature(cf.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor='lightgray'), alpha=.5, zorder=0)
        if polygons is not None:
            polygons.plot(facecolor='none', edgecolor='white', ax=ax)
        if extent is not None:
            ax.set_extent(extent)
        ax.axis('off')

        # scatter plot
        sct = ax.scatter(x, y, c=performance[metric], cmap=cmap, norm=norm, edgecolor='w', lw=1, 
                          s=s, alpha=alpha)
        # # setup: color bar, title...
        cbar = plt.colorbar(sct, ax=ax, shrink=.66)#, orientation='horizontal')
        cbar.set_label(metric, rotation=90)
    
    if 'title' in kwargs:
        fig.text(.5, 1.0, kwargs['title'], ha='center', va='bottom', fontsize=12);
    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight')
        
        
def plot_iterations(
    iters: pd.DataFrame, 
    pareto: pd.DataFrame, 
    best_iter: int, 
    cols: List = ['like1', 'like2'], 
    save: Union[str, Path] = None, 
    **kwargs
):
    """It creates a scatter plot that shows the performance of the iterations in the calibration. On top of the scatter plot a line depicts the Pareto front, from which the best iteration is taken.
    
    Inputs:
    -------
    iters: pd.DataFrame
        A table that contains the two objective functions ("cols") to be plotted
    pareto: pd.DataFrame
        A table that contains the Pareto front from "iters[cols]"
    best_iter: int
        The index of "iters" that contains the best iteration
    save: Union[str, Path]
        If provided, where to save the plot
    """
    
    figsize = kwargs.get('figsize', (4, 4))
    vmax = kwargs.get('vmax', None)
    xlabel = kwargs.get('xlabel', r'$L_{outflow}$')
    ylabel = kwargs.get('ylabel', r'$L_{storage}$')
    title = kwargs.get('title', None)
    
    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(iters[cols[0]], iters[cols[1]], s=1, c='gray', label='iteration', alpha=.2)
    ax.scatter(*iters.loc[best_iter, cols], c='steelblue', label='optimum', s=4)
    ax.plot(pareto[cols[0]], pareto[cols[1]], c='k', lw=1, ls=':', label='pareto front', zorder=0)
    ax.set(xlim=(-.025, vmax),
           xlabel=xlabel,
           ylim=(-.025, vmax),
           ylabel=ylabel)
    if 'title' in kwargs:
        ax.set_title(kwargs['title'])
    fig.legend(frameon=False, loc=1, bbox_to_anchor=[1.175, .7, .1, .2])
    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight');
        
        
def boxplot_parameters(
    parameters: xr.Dataset,
    parameter_range: Optional[Dict[str, List]] = None,
    save: Optional[Union[str, Path]] = None,
    **kwargs
):
    """It creates boxplots comparing the reservoir parameters of a specific model over different runs and reservoirs. 
    One plot is created for each reservoir parameter, in which every box represents the variability of that parameter among reservoirs for a specific run.
    
    Parameters:
    -----------
    parameters: xarray.Dataset
        It contains the reservoir parameters used in different runs and reservoirs. The variables are the reservoir parameters, and it has two dimensions: run and reservoir ID
    parameter_range: dictionary (optional)
        It contains for each parameter (key) the search range during the calibration
    save: string or pathlib.Path (optional)
        If provided, the plot will be saved in this file
        
    Keyword arguments:
    ------------------
    axsize: List
        Size of each of the individual plots
    color: string
        Color of the boxes
    alpha: float
        Transparency of the boxes
    """
    
    axsize = kwargs.get('axsize', (4.5, 4))
    color = kwargs.get('color', 'lightsteelblue')
    alpha = kwargs.get('alpha', 1.0)
    
    ncols = len(parameters)
    fig, axes = plt.subplots(ncols=ncols, figsize=(ncols * axsize[0], axsize[1]), sharey=True)
    if not isinstance(axes, np.ndarray):
        axes = np.array([axes])
    for ax, (parname, da) in zip(axes, parameters.items()):

        # data frame of the parameter
        df = da.to_pandas().transpose()
        df.dropna(inplace=True)

        ax.boxplot(df,
                   # positions=[i + (j - 1) * w],
                   # widths=w * .8,
                   vert=False,
                   patch_artist=True,
                   boxprops=dict(facecolor=color, edgecolor='none', alpha=alpha), 
                   medianprops={'color': 'dimgray'},
                   whiskerprops=dict(color='dimgray', linestyle='-'),
                   showcaps=False,
                   flierprops=dict(marker='.', markersize=2)
                  );
        if parameter_range is not None:
            # ax.set_xlim(parameter_range[parname])
            for v in parameter_range[parname]:
                ax.axvline(v, ls=':', lw=.5, c='k');
            
        ax.set_xlabel(parname)
        ax.tick_params(axis='y', length=0)
        ax.spines[['top', 'right', 'left']].set_visible(False)

    axes[0].set_yticks(range(1, len(df.columns) + 1))
    axes[0].set_yticklabels(df.columns);
    
    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight')
        
        
def compare_attributes(
    df: pd.DataFrame, 
    thr: float, 
    vmin: float = None, 
    vmax: float = None, 
    **kwargs
):
    """Pair plot comparing the attribute values in different data sources
    
    Parameters:
    -----------
    df: pandas.DataFrame
        Values of an attributes in different data sources (columns)
    thr: float
        Minimum value of the attribute used in the selection
    vmin: float
        Minimum attribute value to be plotted
    vmax: float
        Maximum attribute value to be plotted
    """
    
    figsize = kwargs.get('figsize', (4, 4))
    scale = kwargs.get('scale', 'log')
    
    cols = df.columns
    ncols = len(cols) - 1
    
    fig, axes = plt.subplots(ncols=ncols, nrows=ncols, figsize=(ncols * figsize[0], ncols * figsize[1]), sharex=True, sharey=True)   
    if ncols == 1:
        axes = np.array([[axes]])
        
    for i, ax in enumerate(axes.flatten()):
        r = int(i / ncols)
        c = i % ncols
        if c > r:
            ax.axis('off')
            continue
        colx = cols[c]
        coly = cols[r + 1]
        ax.plot([vmin, vmax], [vmin, vmax], c='k', lw=.5, zorder=0)
        ax.vlines(thr, vmin, thr, color='k', ls='--', lw=.5, zorder=0)
        ax.hlines(thr, vmin, thr, color='k', ls='--', lw=.5, zorder=0)
        ax.scatter(df[colx], df[coly], s=10, alpha=.5)
        ax.set_xscale(scale)
        ax.set_yscale(scale)
        if c == 0:
            ax.set_ylabel(coly)
        if r == ncols - 1:
            ax.set_xlabel(colx)

        ax.set(
            xlim=(vmin, vmax),
            ylim=(vmin, vmax),
        );
    
    if 'title' in kwargs:
        fig.suptitle(kwargs['title']);


def boxplot_comparison(
    performance: xr.Dataset,
    ax_dim: str,
    col_dim: str,
    metric: str = 'KGE',
    save: Optional[Union[str, Path]] = None,
    **kwargs,
):
    """
    Generate side-by-side boxplots comparing model performance across different metrics, storage/outflow components,
    and categories (e.g., model types or parameter sets).

    Parameters
    ----------
    performance : xr.Dataset
        An xarray dataset containing performance scores. It must include the specified `ax_dim`, `col_dim`, 
        and a 'metric' dimension, with variables named 'storage' and 'outflow'.
    ax_dim : str
        Dimension of the dataset to map to different subplots (columns of subplots). Typically a grouping like 
        region or threshold.
    col_dim : str
        Dimension that identifies the groups within each boxplot (e.g., different models, scenarios, etc.).
    metric : str, optional
        Performance metric to be plotted (default is 'KGE'). This selects a slice from the 'metric' dimension.
    save : str or Path, optional
        Path to save the figure. If None, the plot is not saved.

    **kwargs : dict, optional
        Additional plot customization options:
        - width (float): Width of each boxplot group (default: 0.15)
        - figsize (tuple): Size of the figure in inches (default: (6, 3))
        - alpha (float): Transparency of the boxplot fill (default: 0.7)

    Notes
    -----
    The function also computes a composite performance index called 'storage & outflow', calculated as:
        1 - sqrt[(1 - storage)^2 + (1 - outflow)^2]
    This is included alongside the individual 'storage' and 'outflow' scores.

    A color-coded legend is added automatically based on the values in `col_dim`.

    Returns
    -------
    None
        Displays the plot and optionally saves it to a file.
    """

    w = kwargs.get('width', .15)
    figsize = kwargs.get('figsize', (6, 3))
    alpha = kwargs.get('alpha', .7)
    ylabel = kwargs.get('ylabel', metric)
    
    colors = ['grey', 'salmon', 'gold', 'steelblue', 'olivedrab']
    colors = {str(key): color for key, color in zip(performance[col_dim].data, colors)}
    
    n = len(performance[ax_dim])
    fig, axes = plt.subplots(ncols=n, figsize=(figsize[0] * n, figsize[1]))#, sharey=True)

    for ax, title in zip(axes, performance[ax_dim].data):
    
        perf = performance.sel({ax_dim: title, 'metric': metric})
        perf = perf.dropna(col_dim, how='all')
        
        perf_dct = {var: perf[var].to_pandas().transpose() for var in ['outflow', 'storage']}
        perf_dct['outflow &\nstorage'] = 1 - ((1 - perf_dct['outflow'])**2 + (1 - perf_dct['storage'])**2)**.5
    
        ticks_labels = {x: var for x, var in enumerate(perf_dct, start=1)}
        for x, var in ticks_labels.items():
            df = perf_dct[var]
            pos = x - (df.shape[1] / 2 - .5) * w
            for i, (col, c) in enumerate(colors.items()):
                if col not in df.columns:
                    continue
                ax.boxplot(
                    df[col],
                    positions=[pos + i * w],
                    widths=[w * .9],
                    patch_artist=True,
                    showfliers=False,
                    capprops=dict(linewidth=0),
                    boxprops=dict(facecolor=c, edgecolor='none', alpha=alpha),
                    whiskerprops=dict(color=c),
                    medianprops=dict(color='k')
                )
        ax.tick_params(axis='x', length=0)
        ax.set(
            ylim=(-1, 1),
            ylabel=ylabel,
            title=title
        )
        ax.spines[['top', 'right', 'bottom']].set_visible(False)
        ax.set_xticks(list(ticks_labels.keys()))
        ax.set_xticklabels(list(ticks_labels.values()))
        ax.set_yticks([-1, -.5, 0, .5, 1])
        
    # Add legend
    legend_handles = [Patch(facecolor=c, edgecolor='none', alpha=.7, label=col) for col, c in colors.items()]
    fig.legend(handles=legend_handles, frameon=False, loc=6, bbox_to_anchor=[.90, 0, .05, 1]);

    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight');


def swarmplot_comparison(
    performance: xr.Dataset,
    ax_dim: str,
    col_dim: str,
    kind: Literal['swarm', 'strip'] = 'swarm',
    metric: str = 'KGE',
    save: Optional[Union[str, Path]] = None,
    **kwargs,
):
    """
    Generate side-by-side plots (swarmplot or stripplot) comparing model performance
    across different metrics, storage/outflow components, and categories (e.g., model types).

    Parameters
    ----------
    performance : xr.Dataset
        An xarray dataset containing performance scores.
        It must include the specified `ax_dim`, `col_dim`, and a 'metric' dimension.
        It should also contain variables named 'storage' and 'outflow'.
    ax_dim : str
        Dimension of the dataset to map to different subplots (columns of subplots).
        Typically a grouping like region or threshold, defining the individual plot columns.
    col_dim : str
        Dimension that identifies the groups within each plot, represented by different
        colors and plotted as individual point swarms/strips (e.g., different models, scenarios).
    kind : {'swarm', 'strip'}, default='swarm'
        Type of plot to generate.
        - 'swarm': Uses `seaborn.swarmplot` to ensure points do not overlap.
        - 'strip': Uses `seaborn.stripplot` which allows jittering to prevent overlap.
    metric : str, default='KGE'
        Performance metric to be plotted (e.g., 'KGE', 'NSE', 'RMSE').
        This selects a specific slice from the 'metric' dimension in the dataset.
    save : str or Path, optional
        Path where the generated figure will be saved. If `None` (default), the plot
        is displayed but not saved to a file.
    **kwargs : dict, optional
        Additional plot customization options that are passed to the function:
        - `figsize` (tuple, default=(20, 3)): Size of the overall figure in inches.
        - `alpha` (float, default=1): Transparency of the swarmplot/stripplot points.
        - `jitter` (float or bool, default=True): Only applicable if `kind='strip'`.
          Determines the amount of jittering applied to points. `True` for default jitter,
          `False` for no jitter, or a float for a specific amount.
        - `linewidth` (float, default=1): Line width for boxplot elements (edges, whiskers, median).
        - `size` (float, default=1.5): Size of the individual points in the swarmplot/stripplot.
        - `width` (float, default=0.5): Width of each box in the boxplot.
        - `width_ratio` (float, default=0.5): Ratio of the width of the empty space between
          `ax_dim` groups compared to the width of a single subplot.
        - `wspace` (float, default=0.25): Horizontal spacing between subplots.
        - `xlim` (tuple, default=(-1, 1)): Tuple specifying the x-axis limits.
        - `ylim` (tuple, default=(-1, 1)): Tuple specifying the y-axis limits.

    Notes
    -----
    The function generates three types of performance visualizations for each `ax_dim` group:
    'outflow', 'storage', and a composite 'outflow & storage'. The composite metric is calculated as:
    $1 - \sqrt{(1 - \text{outflow})^2 + (1 - \text{storage})^2}$

    A color-coded legend is automatically added at the bottom of the figure,
    based on the unique values found in the `col_dim` dimension.

    Returns
    -------
    None
        The function displays the plot and optionally saves it to the specified file path.
    """

    if kind not in ['swarm', 'strip']:
        raise ValueError(f'The attribute "kind" must be either "swarm" or "strip", but {kind} was provided')

    # extract keyword arguments
    figsize = kwargs.get('figsize', (20, 3))
    alpha = kwargs.get('alpha', 1)
    jitter = kwargs.get('jitter', True)
    lw = kwargs.get('linewidth', 1)
    s = kwargs.get('size', 1.5)
    w = kwargs.get('width', .5)
    wratio = kwargs.get('width_ratio', .5)
    wspace = kwargs.get('wspace', 0.25)
    xlim = kwargs.get('xlim', (-2, None))
    ylabel = kwargs.get('ylabel', metric)
    ylim = kwargs.get('ylim', (-1.05, 1.05))
    
    colors = ['grey', 'salmon', 'gold', 'steelblue', 'olivedrab']
    colors = {str(key): color for key, color in zip(performance[col_dim].data, colors)}

    # setup the axes
    n_ax = len(performance[ax_dim])
    n_label = 3
    n_col = 0
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(
        nrows=1, 
        ncols=n_ax * (n_label + 1) - 1, 
        width_ratios=(([1] * n_label + [wratio]) * n_ax)[:-1]
    )
    plt.subplots_adjust(wspace=wspace)

    for i, title in enumerate(performance[ax_dim].data):
    
        perf = performance.sel({ax_dim: title, 'metric': metric})
        
        perf_dct = {var: perf[var].to_pandas().transpose() for var in ['outflow', 'storage']}
        perf_dct['outflow &\nstorage'] = 1 - ((1 - perf_dct['outflow'])**2 + (1 - perf_dct['storage'])**2)**.5
    
        labels = ['outflow', 'storage', 'outflow &\nstorage']
        for j, label in enumerate(labels):

            # create axis
            pos = i * (n_ax + 1) + j
            ax = plt.subplot(gs[pos])

            # boxplot
            box = sns.boxplot(
                perf_dct[label],
                width=w,
                showcaps=False,
                showfliers=False,
                boxprops=dict(facecolor='none', edgecolor='k', linewidth=lw),
                whiskerprops=dict(color='k', linewidth=lw),
                medianprops=dict(color='k', linewidth=lw * 1.5),
                zorder=1,
                ax=ax
            )
                
            # swarmplot
            if kind == 'swarm':
                swarm = sns.swarmplot(
                    perf_dct[label],
                    order=colors.keys(),
                    palette=colors.values(),
                    size=s,
                    alpha=alpha,
                    zorder=0,
                    ax=ax
                )
            elif kind == 'strip':
                strip = sns.stripplot(
                    perf_dct[label],
                    order=colors.keys(),
                    palette=colors.values(),
                    jitter=jitter,
                    size=s,
                    alpha=alpha,
                    zorder=0,
                    ax=ax
                )

            n_col = max(n_col, perf_dct[label].shape[1])
            
            # axis setup
            ax.tick_params(axis='x', length=0)
            ax.set(
                xlim=xlim,
                xlabel=label,
                xticks=[],
                ylim=ylim,                
            )
            ax.spines[['top', 'right', 'bottom']].set_visible(False)
            ax.spines['left'].set_bounds(-1, 1)
            if j % 3 == 0:
                ax.set(
                    ylabel=ylabel,
                    # yticks=[-1, -.5, 0, .5, 1],
                )
                ax.text(-.5, 1.1, f'{chr(97 + i)})', transform=ax.transAxes)
            else:
                ax.set(
                    ylabel=None,
                    yticks=[],
                )
                ax.spines['left'].set_visible(False)
            if j == 1:
                ax.set_title(title)
        
    # Add legend
    legend_handles = [Patch(facecolor=c, edgecolor='none', alpha=.7, label=col) for col, c in colors.items()]
    fig.legend(
        handles=legend_handles, 
        frameon=False, 
        ncol=n_col, 
        loc='lower center', 
        bbox_to_anchor=[.3, -0.175, .4, 0.1]
    )

    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight');


def plot_timeseries(
    ts: pd.DataFrame, 
    cap_mcm: Optional[float] = None, 
    elev_masl: Optional[None] = None,
    dam_hgt_m: Optional[None] = None,
    save: Optional[Union[str, Path]] = None,
    **kwargs
):
    """
    Plots reservoir elevation-storage curve and time series of storage.

    The function creates a two-panel figure:
    - Left panel: elevation vs. storage scatter plot (colored by time).
    - Right panel: time series of storage (also colored by time).

    Optional lines can be drawn to indicate reservoir capacity, dam crest elevation, 
    and base elevation (crest minus dam height).

    Parameters
    ----------
    ts : pandas.DataFrame
        Time series with at least a 'storage' column (in cubic meters).
        If an 'elevation' column is present, it is used in the elevation-storage plot.
    cap_mcm : float, optional
        Reservoir capacity in million cubic meters (hm³). Draws horizontal dashed lines at this value.
    elev_masl : float, optional
        Dam crest elevation in meters above sea level (masl). Draws vertical dashed line at this elevation.
    dam_hgt_m : float, optional
        Dam height in meters. If provided with `elev_masl`, another vertical line is drawn at the base elevation.
    save : str or pathlib.Path, optional
        If provided, path to save the figure as a PNG file. Otherwise, the figure is shown interactively.
    **kwargs : dict, optional
        Additional keyword arguments. Recognized key:
            - 'title': str
                Title for the figure.

    Returns
    -------
    None
        Displays or saves the generated figure.

    Notes
    -----
    - Storage values are converted from cubic meters to million cubic meters (hm³).
    - Points are colored by date using the `coolwarm_r` colormap.
    """

    fig = plt.figure(figsize=(15, 3))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 4], wspace=.05)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1], sharey=ax0)
    
    # reservoir curve
    if 'elevation' in ts.columns:
        ax0.scatter(
            x=ts.elevation, 
            y=ts.storage * 1e-6, 
            c=ts.index,
            cmap='coolwarm_r',
            s=1,
        )
    ax0.set(
        xlabel='elevation (masl)',
        ylabel='storage (hm3)'
    )
    if cap_mcm is not None:
        ax0.axhline(cap_mcm, c='k', ls='--', lw=.5)
        ax1.axhline(cap_mcm, c='k', ls='--', lw=.5)
    if elev_masl is not None:
        ax0.axvline(elev_masl, c='k', ls='--', lw=.5)
        if dam_hgt_m is not None:
            ax0.axvline(elev_masl - dam_hgt_m, c='k', ls='--', lw=.5)
    
    # storage time series
    ax1.scatter(
        x=ts.index,
        y=ts.storage * 1e-6,
        c=ts.index,
        cmap='coolwarm_r',
        s=1
    )
    ax1.set(
        xlabel='date',
        xlim=(ts.first_valid_index(), ts.last_valid_index())
    )
    ax1.yaxis.set_tick_params(labelleft=False)
    # ax1.set_yticklabels([])

    if 'title' in kwargs:
        fig.suptitle(kwargs['title']);

    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight')
        plt.close(fig)


def plot_pareto_front(
    kge_storage: pd.Series, 
    kge_outflow: pd.Series,
    save: Optional[Union[str, Path]] = None,
    **kwargs):
    """
    Plot the Pareto front of KGE scores for storage and outflow.

    This function visualizes the trade-off between two performance metrics 
    (KGE for storage and KGE for outflow) across multiple iterations or 
    parameter sets. Each point is colored by its bivariate KGE score 
    (the Euclidean distance from the ideal point [1,1]). The Pareto front 
    is overlaid, and the best overall iteration (with highest bivariate 
    KGE) is highlighted.

    Parameters
    ----------
    kge_storage : pd.Series
        KGE values for storage, indexed by iteration or model ID.
    kge_outflow : pd.Series
        KGE values for outflow, indexed by iteration or model ID.
    save : str or Path, optional
        Path to save the plot. If None, the plot is shown inline or left open.
    **kwargs :
        Additional keyword arguments, such as:
            - alpha (float): Transparency of scatter points (default is 0.7).
            - figsize (tuple): Figure size in inches (default is (6, 5)).
            - size (float): Size of scatter points (default is 4).
            - title (str): Title for the plot (default is None).

    Returns
    -------
    None
        The function creates a matplotlib plot and optionally saves it to disk.
    """

    alpha = kwargs.get('alpha', .7)
    figsize = kwargs.get('figsize', (6, 5))
    size = kwargs.get('size', 4)
    title = kwargs.get('title', None)
    
    fig, ax = plt.subplots(figsize=figsize)

    # scatter plot of iterations
    kge_2var = 1 - np.sqrt((1 - kge_storage)**2 + (1 - kge_outflow)**2)
    sct = ax.scatter(
        kge_storage, 
        kge_outflow, 
        c=kge_2var, 
        cmap='coolwarm_r', 
        vmin=-1, 
        vmax=1,
        s=size, 
        alpha=alpha
    )
    cbar = plt.colorbar(sct, label='KGE bivariate', shrink=.66)

    # pareto front
    mask_pareto = is_pareto_efficient(kge_storage, kge_outflow)
    pareto_front = pd.concat([kge_storage[mask_pareto], kge_outflow[mask_pareto]], axis=1)
    pareto_front.sort_values(pareto_front.columns[0], inplace=True)
    ax.plot(
        pareto_front.iloc[:, 0], 
        pareto_front.iloc[:, 1], 
        c='k', 
        lw=.8, 
        zorder=2
    )

    # best iteration
    best_iter = kge_2var.idxmax()
    ax.scatter(
        kge_storage[best_iter],
        kge_outflow[best_iter],
        marker='+',
        c='k',
        zorder=3
    )

    # setup
    vlim = (-.025, 1.025)
    ax.plot(vlim, vlim, c='k', ls='--', lw=.5)
    ax.set(
        xlim=vlim,
        xlabel='KGE storage',
        ylim=vlim,
        ylabel='KGE outflow',
        title=title,
    )

    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight')
        plt.close(fig)