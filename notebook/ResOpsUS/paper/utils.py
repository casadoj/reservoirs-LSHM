from typing import Union
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns
import cartopy.crs as ccrs
import pyproj

def map_reservoirs(
    geometry,
    volume: pd.Series,
    area: pd.Series,
    elevation: xr.DataArray,
    uparea: xr.DataArray,
    proj = ccrs.PlateCarree(),
    save: Union[str, Path] = None,
    **kwargs
    ):
    """
    Generate a spatial visualization of reservoirs with elevation and drainage basemaps.

    The function creates a map where reservoir locations are represented by bubbles. 
    The bubble size indicates storage volume, and the bubble color indicates the 
    catchment area. It includes two legends for interpretation.

    Parameters
    ----------
    geometry : GeoDataFrame or similar
        Geospatial data containing the reservoir locations (expects .geometry.x/y).
    volume : pd.Series
        Storage volume data indexed to match the reservoirs. 
        Scaled by sqrt(volume) * 1.5 for marker size.
    area : pd.Series
        Catchment area data indexed to match the reservoirs.
        Scaled by area^0.25 for marker color.
    elevation : xr.DataArray
        Digital Elevation Model (DEM) for the background, plotted in grayscale.
    uparea : xr.DataArray
        Upstream drainage area (accumulation) for highlighting river networks.
    proj : cartopy.crs.Projection or str, optional
        The map projection to use. Can be a Cartopy CRS object or a string 
        identifiable by pyproj. Defaults to ccrs.PlateCarree().
    save : Union[str, Path], optional
        Path where the resulting figure will be saved. If None, the figure 
        is not saved. Defaults to None.
    **kwargs : dict, optional
        Additional configuration for the plot:
        * extent : list of float
            The [lon_min, lon_max, lat_min, lat_max] of the map. 
            Default [-125, -66, 25, 50] (CONUS).
        * figsize : tuple of int
            The width and height of the figure in inches. Default (16, 8).
        * scale : float
            Factor to scale the catchment area colors. Default 1e-3.

    Returns
    -------
    None
        The function renders the plot to the current Matplotlib backend and 
        optionally saves to disk.
    """
    
    cmap = kwargs.get('cmap', 'viridis')
    extent = kwargs.get('extent', None)
    figsize = kwargs.get('figsize', (16, 8))
    scale = kwargs.get('scale', 1e-3)

    # set up plot
    if isinstance(proj, str):
        proj = ccrs.Projection(pyproj.CRS(proj))
    fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': proj})        
    transform = None if proj == ccrs.PlateCarree() else ccrs.PlateCarree()

    # base layers
    kwargs_map = {'ax': ax, 'rasterized': True, 'add_colorbar': False, 'transform': transform}
    elevation.plot(
        cmap='Greys', 
        alpha=1, 
        vmin=-100,
        vmax=4000,
        zorder=0
        )
    uparea.plot(
        cmap='Greys_r', 
        alpha=.3, 
        vmin=-22 * 1e9,
        vmax=22* 1e9,
        zorder=1
        )
        
    # reservoirs
    scatter = ax.scatter(
        geometry.x,
        geometry.y,
        transform=transform,
        c=area.loc[geometry.index]**.25 * scale,
        cmap=cmap,
        s=volume.loc[geometry.index]**.5 * 1.5,
        alpha=.8,
        edgecolor='w',
        lw=.5,
        zorder=10
    )
        
    # color legend
    labels1 = [1000, 10000, 100000, 200000, 400000, 600000]
    handles1, foo = scatter.legend_elements(
        prop='colors', 
        num=list((np.array(labels1))**.25 * scale),
        alpha=.8
    )
    labels1 = [int(label * scale) for label in labels1]
    legend1 = ax.legend(
        handles1,
        labels1,
        title='Catchment\n[10³ km²]', 
        loc='upper left',
        bbox_to_anchor=[1, .5, .1, .35], 
        frameon=False
    )
    fig.add_artist(legend1)

    # size legend
    labels2 = [10, 100, 1000, 10000, 30000]
    handles2, foo = scatter.legend_elements(
        prop='sizes', 
        num=list(np.sqrt(labels2) * 1.5),
        alpha=.8
    )
    legend2 = ax.legend(
        handles2,
        labels2,
        title='Storage\n[hm³]', 
        loc='upper left',
        bbox_to_anchor=[1, .1, .1, .35], 
        frameon=False
    )
    fig.add_artist(legend2)
        
    if extent is not None:
        ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.set_title(None)
    ax.set_aspect('equal')
    ax.axis('off');

    # save
    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight')


def map_performance(
    reservoirs: gpd.GeoDataFrame,
    volume: str,
    performance: str,
    elevation: xr.DataArray,
    uparea: xr.DataArray,
    proj = ccrs.PlateCarree(),
    save: Union[str, Path] = None,
    **kwargs
    ):
    """
    Generate a dual-panel visualization of reservoir performance and spatial distribution.

    The function creates a two-column figure: the left panel shows a combined 
    boxplot and swarmplot of a performance metric, and the right panel shows 
    a spatial map where reservoir markers are sized by volume and colored 
    by performance against elevation and drainage basemaps.

    Parameters
    ----------
    reservoirs : gpd.GeoDataFrame
        Geospatial data containing the reservoir locations, volumes, and 
        performance metrics.
    volume : str
        Column name in `reservoirs` representing storage volume, used to scale 
        marker sizes (scaled by sqrt(volume) * 1.5).
    performance : str
        Column name in `reservoirs` representing the performance metric 
        (e.g., KGE or NSE) for color mapping and statistical plotting.
    elevation : xr.DataArray
        Digital Elevation Model (DEM) for the background, plotted in grayscale.
    uparea : xr.DataArray
        Upstream drainage area (accumulation) for highlighting river networks.
    proj : ccrs.Projection or str, optional
        The map projection to use. Can be a Cartopy CRS object or a string 
        identifiable by pyproj. Defaults to ccrs.PlateCarree().
    save : Union[str, Path], optional
        Path where the resulting figure will be saved. If None, the figure 
        is not saved. Defaults to None.
    **kwargs : dict, optional
        Additional configuration for the plot:
        * alpha : float
            Transparency for the markers and swarmplot. Default 0.8.
        * cmap : str
            Colormap for the performance metric. Default 'coolwarm_r'.
        * extent : list of float
            The [lon_min, lon_max, lat_min, lat_max] of the map. Default None.
        * figsize : tuple of int
            The width and height of the figure in inches. Default (12, 6).
        * lw : float
            Base linewidth for boxes, whiskers, and marker edges. Default 1.2.
        * metric : str
            The label used for the performance metric axis and legend. Default "KGE'".
        * size : float
            Point size for the dots in the swarmplot. Default 4.

    Returns
    -------
    None
        The function renders the plot to the current Matplotlib backend and 
        optionally saves to disk.
    """

    alpha = kwargs.get('alpha', 0.8)
    cmap = kwargs.get('cmap', 'coolwarm_r')
    extent = kwargs.get('extent', None)
    figsize = kwargs.get('figsize', (12, 6))
    lw = kwargs.get('lw', 1.2)
    metric = kwargs.get('metric', "KGE'")
    size = kwargs.get('size', 4)

    # sort GeoDataFrame by decreasing volume
    reservoirs = reservoirs.sort_values(volume, ascending=False)

    # set up plot
    if isinstance(proj, str):
        proj = ccrs.Projection(pyproj.CRS(proj))
    transform = None if proj == ccrs.PlateCarree() else ccrs.PlateCarree()
    fig = plt.figure(figsize=figsize)     
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 6], wspace=.0)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1], projection=proj)

    # swarmplot
    sns.boxplot(
        reservoirs[performance],
        width=.25, 
        showcaps=False,
        showfliers=False,
        boxprops=dict(facecolor='none', edgecolor='k', linewidth=lw),
        whiskerprops=dict(color='k', linewidth=lw),
        medianprops=dict(color='k', linewidth=lw * 1.5),
        zorder=1, 
        ax=ax0
    )
    sns.swarmplot(
        y=reservoirs[performance],
        hue=reservoirs[performance],
        palette=cmap,
        hue_norm=(-1, 1),
        alpha=alpha,
        size=size,
        zorder=0,
        ax=ax0,
        legend=False
    )
    ax0.set(
        ylim=(-1.5, 1.5),
        yticks=[-1, -.5, 0, .5, 1],
        ylabel=metric,
        xticks=[],
    )
    ax0.spines[['top', 'right', 'bottom']].set_visible(False)
    ax0.spines['left'].set_bounds(-1, 1)

    # base layers
    kwargs_map = {'ax': ax1, 'rasterized': True, 'add_colorbar': False, 'transform': ccrs.PlateCarree()}
    elevation.plot(cmap='Greys', vmin=-100, vmax=4000, zorder=0, **kwargs_map)
    uparea.plot(cmap='Greys_r', alpha=.3, vmin=-22 * 1e9, vmax=22* 1e9, zorder=1, **kwargs_map)

    # selected reservoirs
    scatter = ax1.scatter(
        reservoirs.geometry.x,
        reservoirs.geometry.y,
        transform=transform,
        c=reservoirs[performance],
        cmap=cmap,
        vmin=-1,
        vmax=1,
        s=reservoirs[volume]**.5 * 1.5,
        alpha=alpha,
        edgecolor='w',
        lw=lw * .5,
        zorder=2
    )

    # color legend
    cax = fig.add_axes([0.9075, 0.52, 0.008, 0.25])
    bounds = [-1., -.75, -.5, -.25, 0, .25, .5, .75, 1.]
    fig.colorbar(
        scatter,
        cax=cax,
        boundaries=bounds,
        ticks=bounds,
        spacing='proportional',
        shrink=.33
    )
    cax.set_yticks([-1, -.5, 0, .5, 1])
    cax.text(.3, 1.15, metric, ha='left', va='bottom')
    cax.tick_params(left=False)
    for spine in cax.spines.values():
        spine.set_visible(False)

    # size legend
    labels2 = [10, 100, 1000, 10000, 30000]
    handles2, foo = scatter.legend_elements(
        prop='sizes', 
        num=list(np.sqrt(labels2) * 1.5),
        alpha=alpha
    )
    legend2 = ax1.legend(
        handles2,
        labels2,
        title='Storage\n[hm³]',
        loc='center left',
        bbox_to_anchor=[1.02, .3], 
        frameon=False
    )
    fig.add_artist(legend2)
            
    if extent is not None:
        ax1.set_extent(extent, crs=ccrs.PlateCarree())
    ax1.set_title(None)
    ax1.set_aspect('equal')
    ax1.axis('off')

    # save
    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight')