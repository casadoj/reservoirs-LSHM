import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.cm as cm
import matplotlib.colors as colors
from typing import Optional, Union, Dict, List, Tuple
from pathlib import Path

from .inputs import rank_and_filter_data


def plot_release(
    obs_release: pd.Series,
    avg_inflow: float ,
    harmonic: pd.Series,
    linear: xr.DataArray,
    Qmin: Optional[float] = None,
    Qmax: Optional[float] = None,
    save: Optional[Union[str, Path]] = None,
    **kwargs
):
    """
    Plot observed weekly releases alongside harmonic and linear release models.

    This function creates a two-panel plot with the left subplot displaying observed weekly 
    releases as scatter points and the harmonic release as a line plot. The right subplot 
    displays the linear release model results for different standard storage availability (a_st)
    values. Both subplots share the same y-axis. The figure size, title, marker size, marker
    transparency, line width, and colormap can be customized through keyword arguments.

    Parameters:
    obs_release (pd.Series): A Series containing weekly release observations indexed by
        epidemiological week numbers ('epiweek').
    avg_inflow (float): The average inflow (MCM/week).
    harmonic (pd.Series): A Series representing the harmonic standard release data indexed by
        epiweek numbers.
    linear (xr.DataArray): An xarray DataArray representing the linear standard release
        indexed by 'i_st' for weekly standardized inflow and 'a_st' for standard storage availability.
    Qmin (float, optional): The minimum allowable release value for the plot. If provided, a horizontal
        line is added to indicate this value (default: None).
    Qmax (float, optional): The maximum allowable release value for the plot. If provided, a horizontal
        line is added to indicate this value (default: None).
    save (Union[str, Path], optional): The file path or file object to save the figure to. If None,
        the figure is not saved (default: None).

    Keyword Arguments:
    figsize (tuple): The size of the figure in inches. Defaults to (12, 4.5).
    title (str): The title of the plot. If None, no title is set. Defaults to None.
    size (int): The size of the scatter plot markers. Defaults to 8.
    alpha (float): The alpha (transparency) level of the scatter plot markers. Defaults to 0.3.
    linewidth (int): The width of the line in the line plot. Defaults to 2.
    cmap (str): The colormap used for the linear release lines. Defaults to 'Blues'.

    Returns:
    None: The function creates a plot and does not return any value.

    Example Usage:
    >>> plot_release(weekly_obs_series, avg_inflow, harmonic_series, linear_data_array,
                     Qmin=0.5, Qmax=1.5, save="output_plot.png", figsize=(12, 6),
                     title="Weekly Release Comparison", cmap='viridis')

    The function will create a plot with the specified parameters and save it to the path
    provided in the 'save' argument.
    """
    
    figsize = kwargs.get('figsize', (12, 4.5))
    title = kwargs.get('title', None)
    s = kwargs.get('size', 8)
    alpha = kwargs.get('alpha', 0.3)
    lw = kwargs.get('linewidth', 2)
    cmap = kwargs.get('cmap', 'Blues')
    
    fig = plt.figure(figsize=(12, 4.5))
    gs = GridSpec(1, 2, width_ratios=[7.5, 4.5])
    
    # harmonic release
    ax1 = fig.add_subplot(gs[0])
    st_release = (obs_release / avg_inflow) - 1
    ax1.scatter(st_release.index, st_release, c='k', s=s, alpha=alpha, label='observed')
    harmonic.plot(ax=ax1, c='steelblue', lw=lw, label='harmonic')
    if Qmin is not None or Qmax is not None:
        ax1.hlines([Qmin, Qmax], 1, 52, color='steelblue', ls=':', lw=lw/2, label='min-max')
    ax1.legend(loc=2, frameon=False)
    ax1.set(xlim=(1, 52),
           xlabel='epiweek',
           ylabel='standardised release (-)')
    
    # linear release
    ax2 = fig.add_subplot(gs[1], sharey=ax1)
    norm = colors.Normalize(vmin=-.25, vmax=1)
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    for a_st in linear.a_st.data[::5]:
        serie = linear.sel(a_st=a_st).to_pandas()
        color = sm.to_rgba(a_st)
        ax2.plot(serie, color=color, label=r'$A_{{st}}={0:.0f}\%$'.format(a_st * 100))
    ax2.set(xlabel='standardised inflow (-)')
    ax2.legend(loc=2, frameon=False)
    ax2.tick_params(labelleft=None)
    
    if title is not None:
        fig.text(.5, .95, title, fontsize=12, ha='center', va='top')
        
    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight')
        plt.close(fig)
        
        
def plot_nor(
    weekly_storage: pd.DataFrame,
    NOR: pd.DataFrame, 
    n_points: int = 3,
    save: Optional[Union[str, Path]] = None, 
    **kwargs
):
    """
    Plots the weekly storage data with highlighted maximum and minimum observations
    and fills reservoir zones based on NOR (Normal Operating Range) data. The plot
    includes a legend positioned to the right-center of the figure, outside of the axes.
    
    Parameters:
    -----------
    weekly_storage: pandas.DataFrame
        A DataFrame containing the weekly storage data.
            It must contain at least the following columns:
            - epiweek: The epidemiological week of the observation.
            - s_st: The storage percentage for the corresponding epiweek.
    NOR: pandas.DataFrame
        A DataFrame containing the Normal Operating Range data.
            It must contain at least the following columns:
            - index: The index, which represents the epiweek.
            - flood: The flood threshold percentage for the corresponding epiweek.
            - conservation: The conservation threshold percentage for the corresponding epiweek.
    n_points: integer (optional)
        Number of maximum/minimum weekly storage values used to fit the flood/conservative storage harmonic function
    save: string or pathlib.Path (optional)
        Path to save the figure. If None, the figure is not saved.
    
    Keyword Arguments:
    ------------------
    figsize (tuple of int): The size of the figure in inches (default: (6, 4.5)).
    title (str): The title of the plot. If None, no title is set (default: None).
    size (int): The size of the scatter plot markers for max. and min. observations (default: 8).
    alpha (float): The alpha (transparency) level of the scatter plot markers for observations
                   (default: 0.3).
    
    Returns:
    --------
    None: The function creates a plot and does not return any value.
    
    Example Usage:
    --------------
    >>> plot_nor(weekly_storage_df, NOR_df, save="output_plot.png", figsize=(10, 7), title="Weekly Storage vs. NOR")
    
    The function will plot the data, show the legend on the right side, and optionally save the
    figure to the specified path.
    """
    
    figsize = kwargs.get('figsize', (6, 4.5))
    title = kwargs.get('title', None)
    s = kwargs.get('size', 8)
    alpha = kwargs.get('alpha', 0.3)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # observations
    ax.scatter(weekly_storage.epiweek, weekly_storage.s_st, c='k', s=8, alpha=alpha, label='observed')
    top = rank_and_filter_data(weekly_storage, 's_st', n_points, ascending=False)
    ax.scatter(top.epiweek, top.s_st, c='green', s=s, alpha=alpha * 2, label='max. obs.')
    bottom = rank_and_filter_data(weekly_storage, 's_st', n_points, ascending=True)
    ax.scatter(bottom.epiweek, bottom.s_st, c='maroon', s=s, alpha=alpha * 2, label='min. obs.')
    
    # reservoir zones
    ax.fill_between(NOR.index, NOR.flood, 1, color='whitesmoke', zorder=0)
    ax.fill_between(NOR.index, NOR.conservation, NOR.flood, color='lightsteelblue', zorder=0, label='NOR')
    ax.fill_between(NOR.index, 0, NOR.conservation, color='whitesmoke', zorder=0)
    
    if title is not None:
        ax.set_title(title)
    ax.text(0.01, 0.975, 'flood pool', va='top', transform=ax.transAxes)
    ax.text(0.01, 0.025, 'conservation pool', va='bottom', transform=ax.transAxes)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
    ax.set(xlim=(1, 52),
           xlabel='epiweek',
           ylim=(0, 1),
           ylabel='reservoir filling (-)');
    
    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight')
        plt.close(fig)
        
        
def find_closest_dam():
    """Finds the dam that is closest in terms of purposes served and Euclidean distance
    
    Parameters:
    -----------
    dam_attr attributes of target dam
    other_dams table of attributes for possible canditate dams to replicate
    distance_only allows for use of closest distance only (disregarding purpose)
    
    Returns:
    --------
    GRAND_ID of the target dam
    """
    
    pass


def aggregate_to_epiweeks(daily: pd.DataFrame) -> pd.DataFrame:
    """Aggregate daily timeseries to epiweek
    
    Parameters:
    -----------
    daily: pandas.DataFrame
        Daily time series that includes fields 'epiweek', 'i' inflow, 's' storage, 'r' release
        
    Returns:
    --------
    weekly: pandas.DataFrame
        Weekly time series
    """
    
    # find the first timestep that represents the beginning of a week
    start_snip = daily.index.get_loc(daily[(daily['epiweek'].diff() != 0)].index[0])

    # Check if the first water week duration is greater than 7 days
    if start_snip not in range(1, 8):
        raise ValueError("first water week duration > 7 days!")

    # snip the data if necessary
    if start_snip < 7:
        daily_snipped = daily.iloc[start_snip:].copy()
    else:
        daily_snipped = daily.copy()

    # Perform aggregation
    daily_snipped['s_end'] = daily_snipped['s'].shift(-7)
    weekly = daily_snipped.groupby(['year', 'epiweek']).agg(
        i=('i', 'sum'),
        r=('r', 'sum'),
        s_start=('s', 'first'),
        s_end=('s_end', 'first')
    ).reset_index()

    # Filter out epiweek 53 and rows where s_end is NA
    weekly = weekly[(weekly['epiweek'] != 53) & (~weekly['s_end'].isna())]

    return weekly


def epiweek_to_date(year, epiweek):
    # Convert epiweek to a date string where the week starts on Sunday
    d = f'{year:.0f}-W{int(epiweek)-1}-0'  # Set to Sunday
    # Use '%U' for weeks starting on Sunday and '%w' set to '0' for Sunday
    return pd.to_datetime(d, format='%Y-W%U-%w')


def back_calc_missing_flows(
    weekly: pd.DataFrame,
    min_weeks: int = 260,
) -> pd.DataFrame:
    """Compute inflow or release from mass balance (if either is missing)
    
    Parameters:
    -----------
    weekly: pandas.DataFrame
        Weekly time series obtained from function `aggregate_to_epiweeks()`
    min_weeks: integer
        Minimum allowable data points (weeks) to use release and inflow without any back-calculating. Default value set to 5 years
    
    Returns:
    --------
    weekly: pandas.DataFrame
        Weekly time series filled in by mass balance in the 'min_weeks' condition is not met
    """

    # Compute the change in storage and back-calculate release (r_) and inflow (i_)
    weekly['s_change'] = weekly['s_end'] - weekly['s_start']
    weekly['r_'] = np.where((weekly['i'] - weekly['s_change']) < 0, 0, weekly['i'] - weekly['s_change'])
    weekly['i_'] = weekly['r'] + weekly['s_change']

    # Filter to full data points where neither r nor i is NA
    full_data_points = weekly.dropna(subset=['r', 'i'])

    # Check the number of data points on the most data-scarce epiweek
    if full_data_points.shape[0] == 0:
        data_points_on_most_data_scarce_epiweek = -np.inf
    else:
        data_points_on_most_data_scarce_epiweek = full_data_points.groupby('epiweek').size().min()

    # back-calculate if there aren't enough data points
    if data_points_on_most_data_scarce_epiweek < min_weeks:
        
        # Count missing values
        missing_i = weekly['i'].isna().sum()
        missing_r = weekly['r'].isna().sum()

        # Decide which variable to back-calculate based on which has fewer missing values
        if missing_i <= missing_r:
            weekly['i'] = np.where(weekly['i'].isna() & ~weekly['r'].isna(), weekly['i_'], weekly['i'])
            weekly['r'] = np.where(weekly['r_'].isna(), weekly['r'], weekly['r_'])
        else:
            weekly['r'] = np.where(weekly['r'].isna() & ~weekly['i'].isna(), weekly['r_'], weekly['r'])
            weekly['i'] = np.where(weekly['i_'].isna(), weekly['i'], weekly['i_'])

    return weekly[['year', 'epiweek', 's_start', 'i', 'r']]