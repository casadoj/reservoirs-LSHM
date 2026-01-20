import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from timezonefinder import TimezoneFinder
from pytz.exceptions import AmbiguousTimeError, NonExistentTimeError
from typing import Literal, Optional


def convert_to_utc(
    lon: float, 
    lat: float, 
    series: pd.DataFrame
) -> pd.DataFrame:
    """Takes a time series with no timezone definition, finds the timezone based on the
    coordinates of the location, and converts the datetime index to UTC

    Parameters:
    -----------
    lon: float
        The longitude of the location for which the timezone is to be determined.
    lat: float
        The latitude of the location for which the timezone is to be determined.
    series: pandas.DataFrame
        The input time series. Its index is expected to be a pandas.DatetimeIndex.

    Returns:
    --------
    pandas.DataFrame: 
        A new pandas.DataFrame with the same data but with its DatetimeIndex converted to UTC.
    """

    # find timezone
    tf = TimezoneFinder()
    local_timezone = tf.timezone_at(lng=lon, lat=lat)
    

    # define timezone
    ts_local = series.copy()
    if ts_local.index.tz is None:
        try:
            ts_local.index = ts_local.index.tz_localize(local_timezone)
        except NonExistentTimeError as e:
            # Handle DST gaps (e.g., clocks "spring forward")
            ts_local.index = ts_local.index.tz_localize(local_timezone, nonexistent='shift_forward')
        except AmbiguousTimeError as e:
            # Handle ambiguous times (e.g., clocks "fall back")
            ts_local.index = ts_local.index.tz_localize(local_timezone, ambiguous=False)
    else:
        ts_local.index = ts_local.index.tz_convert(local_timezone)

    # convert to UTC with offset
    ts_utc = ts_local.copy()
    ts_utc.index = ts_utc.index.tz_convert('UTC')
    
    # store timezone and offset for optional use
    convert_to_utc.local_tz = local_timezone
    convert_to_utc.offset = ts_utc.index.hour.min()

    return ts_utc


def water_balance_iterative(
    Vo: float, 
    inflow: pd.Series, 
    outflow: pd.Series
) -> pd.Series:
    """Calculates by water balance the time series of reservoir storage given the initial storage and the time
    series of inflow and outflow.

    Parameters:
    -----------
    Vo: float
        Initial reservoir storage
    inflow: pandas.Series
        Time series of reservoir inflow. It can't contain missing values
    outflow: pandas.Series
        Time series of reservoir outflow. It can't contain missing values

    Returns:
    --------
    pandas.DataFrame
        Time series of reservoir storage at the beginning of each time step
    """

    # concatenate inputs
    df = pd.concat([inflow, outflow], axis=1)
    df.columns = ['inflow', 'outflow']
    if df.isnull().any().any():
        raise ValueError("The input time series can't contain missing values")
    
    # time step duration
    duration = pd.Series(df.index.diff().total_seconds(), index=df.index).shift(-1)
    water_balance_iterative.duration = duration

    # compute storage by water balance
    V = Vo
    storage = [V]
    for T, I, O in zip(duration, df.inflow, df.outflow):
        V += (I - O) * T
        storage.append(V)
    storage = pd.Series(storage[:-1], index=df.index, name='storage_wb')
        
    return storage


def water_balance_timestep(
    storage: pd.Series,
    inflow: pd.Series, 
    outflow: pd.Series
) -> pd.Series:
    """Calculates by water balance the time series of reservoir storage. It operates independently for each
    time step, using the time series of storage, inflow and outflow. It works if there are missing values in
    the input time series.

    The objective of this function it to check if the convertion of time series into 00 UTC conserves mass.

    Parameters:
    -----------
    storage: pd.Series
        Time series of reservoir storage
    inflow: pandas.Series
        Time series of reservoir inflow
    outflow: pandas.Series
        Time series of reservoir outflow

    Returns:
    --------
    pandas.DataFrame
        Time series of reservoir storage at the beginning of each time step
    """

    # concatenate inputs
    df = pd.concat([inflow, outflow, storage], axis=1)
    df.columns = ['inflow', 'outflow', 'storage']
    
    # empty result object
    result = pd.Series(index=df.index, name='storage_wb')
    
    # time step duration
    duration = pd.Series(df.index.diff().total_seconds(), index=df.index).shift(-1)
    water_balance_timestep.duration = duration
    
    # compute water balance
    mask = df[['inflow', 'outflow', 'storage']].notna().all(axis=1) # time steps with all data
    storage_wb = df.loc[mask, 'storage'] + (df.loc[mask, 'inflow'] - df.loc[mask, 'outflow']) * duration[mask]
    storage_wb = storage_wb.shift().dropna()
    
    # assign water balance to the results
    idx = result.index.intersection(storage_wb.index)
    result[idx] = storage_wb[idx]
    
    # fill in gaps with original data
    missing = df.storage.notnull() & result.isnull()
    result[missing] = df.loc[missing, 'storage']
    
    return result



def reindex_to_00utc(series_utc: pd.DataFrame) -> pd.DataFrame:
    """Takes a daily time series in the UTC timezone, but with an offset (at hours different that 00 am),
    and estimates the associated value at 00 am UTC.

    To shift the time, it uses pandas.reindex() with linear interpolation.
    
    Parameters:
    -----------
    series_utc: pandas.DataFrame
        Daily time series with a DatetimeIndex in UTC

    Returns:
    --------
    pandas.DataFrame
        A new pandas.DataFrame with the same data but shifted to 00 am UTC.
    """
    
    # timestamps with offset UTC
    idx_offset = series_utc.index
    
    # timestamps at 00 UTC
    dates = series_utc.index.date
    offset = series_utc.index.hour.min()
    if offset == 0:
        start, end = series_utc.index[0], series_utc.index[-1]
    if offset < 12:
        start = series_utc.index.min().floor('D')
        end = series_utc.index.max().floor('D')
    elif offset >= 12:
        start = series_utc.index.min().ceil('D')
        end = series_utc.index.max().ceil('D')
    idx_00 = pd.date_range(
        start=start,
        end=end,
        freq='D',
        tz='UTC',
        name=series_utc.index.name
        )
    
    # interpolate linearly the values at 00 UTC
    idx = idx_00.union(idx_offset)
    series_reindex = series_utc.reindex(idx).interpolate(
        'time', 
        limit=1, 
        limit_direction='both',
    )
    reindex_to_00utc.series_reindex = series_reindex
    series_00utc = series_reindex.loc[idx_00]

    # # correct storage with water balance
    # storage_wb = water_balance_timestep(series_00utc.storage, series_00utc.inflow, series_00utc.outflow)
    # mask_error = ~np.isclose(series_00utc.storage, storage_wb)
    # series_00utc.loc[mask_error, 'storage'] = storage_wb[mask_error]

    return series_00utc

    
def resample_to_00utc(
    series_utc: pd.DataFrame, 
    kind: Literal['average', 'instant'] = 'average'
) -> pd.DataFrame:
    """Takes a daily time series in the UTC timezone, but with an offset (at hours different that 00 am),
    and estimates the associated value at 00 am UTC.

    To shift the time, it resamples the data to hourly values and then aggregates them back to daily values
    at 00 am. The procedure differs whether the variable represents average or instantaneous values. 
        'average': the first case, the hourly resampling uses forward filling and the aggregation back to 
                   daily values calculates the mean
        `instant`: the hourly resampling uses linear interpolation and the aggregation back to daily values
                   picks the values at 00 am
    
    Parameters:
    -----------
    series_utc: pandas.DataFrame
        Daily time series with a DatetimeIndex in UTC

    Returns:
    --------
    pandas.DataFrame
        A new pandas.DataFrame with the same data but shifted to 00 am UTC.
    """

    series_utc = series_utc.copy()
    
    if kind not in ['instant', 'average']:
        raise ValueError('"kind" must be one of the following options: "average", "instant"')

    # define extent of the resulting time series
    offset = series_utc.index.hour.min()
    print(f'offset: {offset} h')
    if offset == 0:
        start, end = series_utc.index[0], series_utc.index[-1]
    if offset < 12:
        start = series_utc.index.min().floor('D')
        end = series_utc.index.max().floor('D')
        series_utc.loc[start] = np.nan
    elif offset >= 12:
        start = series_utc.index.min().ceil('D')
        end = series_utc.index.max().ceil('D')
        series_utc.loc[end] = np.nan
    series_utc.sort_index(inplace=True)
        
    # if offset == 0:
    #     pass
    # elif offset < 12:
    #     series_utc.loc[series_utc.index[0].floor('D')] = np.nan
    # elif offset >= 12:
    #     series_utc.loc[series_utc.index[-1].ceil('D')] = np.nan
    # series_utc.sort_index(inplace=True)
    
    # resample hourly values
    if kind == 'instant':
        series_hourly = series_utc.resample('h').interpolate(method='linear', limit=24)
    elif kind == 'average':
        series_hourly = series_utc.resample('h').ffill(limit=24)        
    resample_to_00utc.series_hourly = series_hourly
    
    # daily values
    groupby = series_hourly.groupby(series_hourly.index.date)
    if kind == 'instant':
        # extract instant values at 00 UTC
        series_00 = groupby.first()
    elif kind == 'average':
        # daily mean
        series_00 = groupby.mean()
        
    # # remove values if less than half a day is covered
    # mask = groupby.count() < 12
    # series_00[mask] = np.nan
    
    # defice UTC as the timezone
    series_00.index = pd.DatetimeIndex(series_00.index, tz='UTC')
        
    return series_00.loc[start:end]


def plot_tz_conversion(
    original: pd.DataFrame,
    shifted: pd.DataFrame,
    plot_dates: bool = True,
    **kwargs
):
    """
    Plot and compare original and timezone-shifted time series data 
    for 'inflow', 'outflow', and 'storage' variables.

    This function visualizes the effects of timezone conversion by 
    plotting two DataFrames (one in local time and one in UTC) on 
    separate subplots. Step plots are used for 'inflow' and 'outflow',
    while line plots with markers are used for 'storage'. Vertical lines
    mark day boundaries in the UTC-shifted data.

    Parameters
    ----------
    original : pd.DataFrame
        Time series data in local timezone with columns: 'inflow', 'outflow', 'storage'.
    
    shifted : pd.DataFrame
        Time series data converted to UTC, with the same structure as `original`.
    plot_dates: boolean
        Whether to plot vertical lines to show the dates in UTC.

    Returns
    -------
    None
        Displays a matplotlib figure with three subplots and a shared legend.
    """
    
    original = original.copy()
    original.loc[original.index[-1] + pd.Timedelta(days=1)] = np.nan
    shifted = shifted.copy()
    shifted.loc[shifted.index[-1] + pd.Timedelta(days=1)] = np.nan
    
    fig, axes = plt.subplots(figsize=(12, 9), nrows=3, sharex=True)
    for ax, variable in zip(axes, ['inflow', 'outflow', 'storage']):
        drawstyle = 'steps-post' if variable in ['inflow', 'outflow'] else 'default'
        marker = None if variable in ['inflow', 'outflow'] else '+'
        # original values 
        original[variable].plot(
            ax=ax, 
            label='local timezone', 
            drawstyle=drawstyle,
            marker=marker
        )
        # shifted values
        shifted[variable].plot(
            ax=ax, 
            label='UTC timezone', 
            drawstyle=drawstyle,
            linestyle='--',
            marker=marker
        )
        if plot_dates:
            for date in np.unique(shifted.index.date):
                ax.axvline(date, c='k', lw=.5, ls=':')
        ax.set_title(variable)

    if 'xlim' in kwargs:
        ax.set_xlim(kwargs['xlim'])
    
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, loc=6, bbox_to_anchor=[.91, .4, .1, .3]);