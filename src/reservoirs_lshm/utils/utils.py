import os
os.environ['USE_PYGEOS'] = '0'
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import geopandas as gpd
import cartopy.feature as cfeature
import cartopy.crs as ccrs
from typing import Union, List, Dict, Tuple, Optional, Literal
from pathlib import Path
import xml.etree.ElementTree as ET
from scipy.stats import gumbel_r, gaussian_kde
from tqdm.auto import tqdm



def find_connections(
    dst: gpd.GeoDataFrame,
    src: gpd.GeoDataFrame,
    max_distance: float = 0.01
) -> Dict:
    """Finds a mapping between the indices in two GeoDataFrames based on geographical proximity
    
    Parameters:
    -----------
    dst: geopandas.GeoDataFrame
        Points for which an ID is searched
    src: geopandas.GeoDataFrame
        Points that are used as the source of the ID
        
    Returns:
    --------
    mapping: dictionary
        The keys are indices in 'src', and the values indices in 'dst'
    """
    
    mapping = {}
    for ID, row in tqdm(dst.iterrows(), total=dst.shape[0]):
        # compute "distance" from all points in the source
        diff = ((src.geometry.x - row.geometry.x)**2 + (src.geometry.y - row.geometry.y)**2)**.5
        if diff.min() <= max_distance:
            mapping[ID] = diff.idxmin()          

    return mapping



def filter_domain(
    points: gpd.GeoDataFrame,
    domain: xr.DataArray
):
    """Selects the points inside the model domain.
    
    Parameters:
    -----------
    points: geopandas.GeoDataFrame
    domain: xarray.DataArray
        Boolean map that defines the model domain: 1 defines the domain
    """
    
    def nearest_domain_value(point, domain):
        x_idx = domain.x.get_indexer([point.x], method='nearest')[0]
        y_idx = domain.y.get_indexer([point.y], method='nearest')[0]
        return dom.isel(x=x_idx, y=y_idx).item()
    
    # filter by extent
    lon_min, lat_min, lon_max, lat_max = np.round(domain.rio.bounds(), 6)
    mask_extent = (lon_min <= points.geometry.x) &  (points.geometry.x <= lon_max) & (lat_min <= points.geometry.y) & (points.geometry.y <= lat_max)
    points = points[mask_extent]
    
    # keep points inside the domain
    pbar = tqdm(points.iterrows(), total=points.shape[0])
    mask_domain = [
        ID for ID, point in pbar
        if domain.sel(x=point.geometry.x, y=point.geometry.y, method='nearest').item() == 1
    ]
    points = points.loc[mask_domain]
    
    return points



def filter_reservoirs(
    catchment: pd.Series,
    volume: pd.Series,
    catch_thr: Optional[float] = 10,
    vol_thr: Optional[float] = 10
) -> pd.Series:
    """
    Filters reservoirs based on minimum catchment area and volume thresholds.
    
    Parameters:
    -----------
    catchment: pandas.Series
        Reservoir catchment area
    volume: pandas.series
        Reservoir volume
    catch_thr: float or None
        Minimum catchment area that will be selected. Make sure that the units are the same as "catchment". If "catchment" does not report a value for a reservoir, that reservoir WILL NOT be removed, as this value can be estimated later on
    vol_thr: float or None
        Minimum reservoir volume required for a reservoir to be selected. Make sure that the units are the same as "volume". If "volume" does not report a value for a reservoir, the reservoir WILL be removed
        
    Returns:
    -------
    pd.Series
        A boolean pandas Series where True indicates that a reservoir meets both the catchment and volume thresholds.
    """
    
    assert catchment.shape == volume.shape, '"catchment" and "volume" must have equal shape'
    
    n_reservoirs = catchment.shape[0]
    
    if catch_thr is not None:
        mask_catch = (catchment.isnull()) | (catchment >= catch_thr)
        print('{0} out of {1} reservoirs exceed the minimum catchment area of {2} km2 ({3} missing values)'.format(
            mask_catch.sum(),
            n_reservoirs,
            catch_thr,
            catchment.isnull().sum()
        ))
    else:
        mask_catch = pd.Series(True, index=catchmen.index)
    
    if vol_thr is not None:
        mask_vol = volume >= vol_thr
        print('{0} out of {1} reservoirs exceed the minimum reservoir volume of {2} hm3 ({3} missing values)'.format(
            mask_vol.sum(),
            n_reservoirs,
            vol_thr,
            volume.isnull().sum()
        ))
    else:
        mask_vol = pd.Series(True, index=volume.index)
        
    print('{0} out of {1} reservoirs exceed the minimum catchment area ({2} km2) and the minimum reservoir volume ({3} hm3)'.format(
        (mask_catch & mask_vol).sum(),
        n_reservoirs,
        catch_thr,
        vol_thr
    ))
    
    return mask_catch & mask_vol



def remove_duplicates(
    df: pd.DataFrame,
    duplicates_col: str,
    select_col: str,
    ascending: bool = False,
) -> pd.DataFrame:
    """Given a DataFrame, it identifies duplicate entries in a column and selects that with the largest value in another column

    Parameters:
    -----------
    df: pd.DataFrame
        table from which duplicates will be removed
    duplicates_col: string
        column in "df" where duplicated values will be identified
    select_col: string
        column in "df" used to select one entry from the duplicates. For each duplicated value in "duplicated_col", the largest (ascending=False) or smallest (ascending=True) value in "select_col" will be kept
    ascending: boolean
        whether to sort the 'select_col' in ascending (True) or descending (False) order
        
    Returns:
    --------
    A DataFrame similar to the input, but with duplicates removed
    """
    
    df_ = df.copy()
    
    for value, count in df_[duplicates_col].value_counts().items():
        if count > 1:
            remove_idx = df_.loc[df_[duplicates_col] == value].sort_values(select_col, ascending=ascending).index[1:]
            df_.drop(remove_idx, axis=0, inplace=True)
        else:
            break
            
    return df_
            
            
            
def select_reservoirs(
    df: gpd.GeoDataFrame,
    sort: str,
    storage: str,
    target: float,
    plot: bool = True,
    **kwargs
) -> gpd.GeoDataFrame:
    """Selects reservoirs that fulfil a target total storage capacity by prioritizing based on another characteristic
    
    Inputs:
    -------
    df:    geopandas.GeoDataFrame
        Table of reservoirs
    sort:  string
        Name of the field in 'df' that will be use to sort (prioritize) the selection
    storage: string
        Name of the field in 'df' that contains the reservoir storage capacity
    target: float
        Total storage to be reached
    plot:    boolean
        If True, a map of the selected reservoirs will be plotted. The size of the dots represents the reservoir storage capacity and the colours the sorting field.
    
    Outputs:
    --------
    df_sel: geopandas.DataFrame
        A subset of 'df' with the selection of reservoirs.
    """
    
    mask = df.sort_values(sort, ascending=False)[storage].cumsum() <= target
    df_sel = df.loc[mask]
    volume = df_sel[storage].sum()
    
    if plot:
        fig, ax = plt.subplots(figsize=kwargs.get('figsize', (20, 5)), subplot_kw=dict(projection=ccrs.PlateCarree()))
        ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', edgecolor='face', facecolor='lightgray'), alpha=.5, zorder=0)
        if 'c' in kwargs:
            if isinstance(kwargs['c'], str):
                c = kwargs['c']
            elif isinstance(kwargs['c'], pd.Series):
                c = kwargs['c'].loc[mask]
        else:
            c = df_sel[sort]
        scatter = ax.scatter(df_sel.geometry.x, df_sel.geometry.y, s=df_sel[storage] / 1000, cmap=kwargs.get('cmap', 'coolwarm'), c=c, alpha=kwargs.get('alpha', .5))
        if 'title' in kwargs:
            ax.text(.5, 1.07, kwargs['title'], horizontalalignment='center', verticalalignment='bottom', transform=ax.transAxes, fontsize=12)
        text = '{0} reservoirs   {1:.0f} km³'.format(mask.sum(), volume / 1000)
        ax.text(.5, 1.02, text, horizontalalignment='center', verticalalignment='bottom', transform=ax.transAxes)
        ax.axis('off');
        # if 'c' in kwargs:
        #     if isinstance(kwargs['c'], pd.Series):
        #         legend1 = ax.legend(*scatter.legend_elements(prop='colors', num=4, alpha=.5), title=kwargs.get('legend_title', ''), bbox_to_anchor=[1.025, .65, .09, .25], frameon=False)
        #         ax.add_artist(legend1)
        legend2 = ax.legend(*scatter.legend_elements(prop='sizes', num=4, alpha=.5), title='storage (km³)', bbox_to_anchor=[1.025, .35, .1, .25], frameon=False)
        ax.add_artist(legend2);
    
    return df_sel



def xml_parameters(xml: Union[str, Path], pars: Union[str, List[str]] = None) -> Dict:
    """It extracts the temporal information from the settings XML file.
    
    Input:
    ------
    xml:         Union[str, Path] 
        A XML settings file (path, filename and extension)
    pars:        Union[str, List[str]]
        Name(s) of the parameters to be extracted
        
    Output:
    -------
    parameters:  Dict
        Keys are parameter names and values the calibrated parameter value
    """
    
    # extract temporal info from the XML
    tree = ET.parse(xml)
    root = tree.getroot()
    
    if pars is None:
        pars = ['b_Xinanjiang', 'UpperZoneTimeConstant', 'LowerZoneTimeConstant', 'LZThreshold',
                'GwPercValue', 'GwLoss', 'PowerPrefFlow', 'SnowMeltCoef',
                'AvWaterRateThreshold' , 'LakeMultiplier', 'adjust_Normal_Flood', 'ReservoirRnormqMult', 
                'QSplitMult', 'CalChanMan', 'CalChanMan2', 'ChanBottomWMult', 'ChanDepthTMult', 'ChanSMult']
    
    parameters = {par: float(root.find(f'.//textvar[@name="{par}"]').attrib['value']) for par in pars}
        
    return parameters



def CDF(series: pd.Series):
    """It estimates the value associated to a specific return period based on the observed time series and the Gumbel distribution
    
    Input:
    ------
    series: pd.Series
        Time series from which the annual maxima (therefore the index must be a timestamp) will be extracted and then used to fit a Gumbel distribution
        
    Ouput:
    ------
    CDF: pd.Series
        A series in which the index is the sorted annual maxima and the values the probability of non exceeding that value
    """
    
    # annual maxima
    maxima = series.groupby(series.index.year).max()
    maxima.sort_values(ascending=True, inplace=True)
    
    # fit gumbel distribution
    pars = gumbel_r.fit(maxima.values)
    
    CDF = pd.Series(gumbel_r.cdf(maxima, *pars), index=maxima)
    
    return CDF



def get_normal_value(series: pd.Series):
    """Given values of a variable, it estimates the Gaussian kernel density and ouputs the value of the variable with the highest density.
    
    Input:
    ------
    series: pd.Series
        Values of any variable
        
    Ouput:
    ------
    x: float
        Value of the input variable with the highest Gaussian density.
    """
    
    series_ = series.dropna()
    kde = gaussian_kde(series_)
    x = np.linspace(series_.min(), series_.max(), 1000) #serie.copy().sort_values()
    y = kde(x)
    return x[np.argmax(y)]



def return_period(series: pd.Series, T: float = 100) -> float:
    """It estimates the value associated to a specific return period based on the observed time series and the Gumbel distribution
    
    Input:
    ------
    series: pd.Series
        Time series from which the annual maxima (therefore the index must be a timestamp) will be extracted and then used to fit a Gumbel distribution
    T: int
        Return period (in years) to be estimated
        
    Output:
    -------
    x: float
        Value of the input variable associated with a return period of 'T' years.
    """
    
    series_ = series.dropna()
    
    # annual maxima
    maxima = series_.groupby(series_.index.year).max()
    maxima.sort_values(ascending=True, inplace=True)
    
    # fit gumbel distribution
    pars = gumbel_r.fit(maxima.values)
    return_period.parameters = pars
    
    # discharge associated to return period
    x = gumbel_r.ppf(1 - 1 / T, *pars)
    
    return np.float16(x)



def dict2da(dictionary: Dict, dim: str) -> xr.DataArray:
    """It converts a dictionary of xarray.Datarray into a single xarray.DataArray combining the keys in the dictionary in a new dimension
    
    Inputs:
    -------
    dictionary: dict. A dictionary of xarray.DataArray
    dim:        str. Name of the new dimension in which the keys of 'dictionary' will be combined
    
    Output:
    -------
    array:      xr.DataArray.
    """
    
    if isinstance(dictionary, dict) is False:
        return 'ERROR. The input data must be a Python dictionary.'
        
    data = list(dictionary.values())
    coord = xr.DataArray(list(dictionary), dims=dim)

    return xr.concat(data, dim=coord)



def read_static_map(path: Union[Path, str],
                    x_dim: str = 'lon',
                    y_dim: str = 'lat',
                    crs: str = 'epsg:4326',
                    var: str = 'Band1') -> xr.DataArray:
    """It reads the NetCDF of a LISFLOOD static map as an xarray.DataArray.

    Parameters:
    -----------
    path: Path
        name of the NetCDF file to be opened
    x_dim: str
        name of the dimension that represents coordinate X
    y_dim: str
        name of the dimension that represents coordinate Y
    crs: str
        EPSG code of the coordinate reference system (for instance 'epsg:4326')
    var: str
        name of the variable to be loaded from the NetCDF file

    Returns:
    --------
    xr.DataArray
        a map with coordinates "x_dim", "y_dim" in the reference system "crs"
    """
    
    # load dataset
    da = xr.open_mfdataset(path, chunks=None)[var].compute()

    # set spatial dimensions
    da = da.rio.set_spatial_dims(x_dim=x_dim, y_dim=y_dim)
    
    # define coordinate system
    da = da.rio.write_crs(crs)

    return da


def upstream_pixel(
    lat: float,
    lon: float,
    upArea: xr.DataArray
) -> tuple:
    """This function finds the upstream coordinates of a given point
    
    Parameteres:
    ------------
    lat: float
        latitude of the input point
    lon: float
        longitued of the input point
    upArea: xarray.DataArray
        map of upstream area
        
    Returns:
    --------
    lat: float
        latitude of the inmediate upstream pixel
    lon: float
        longitued of the inmediate upstream pixel
    """
    

    # Determine coordinate system
    lat_coord = 'lat' if 'lat' in upArea.coords else 'y'
    lon_coord = 'lon' if 'lon' in upArea.coords else 'x'
    
    try:

        # upstream area of the input coordinates
        area = upArea.sel({lat_coord: lat, lon_coord: lon}, method='nearest').item()
        
        # spatial resolution of the input map
        resolution = np.mean(np.diff(upArea[lon_coord].values))
        
        # Define window around the input pixel
        window = 1.5 * resolution
        upArea_window = upArea.sel({lat_coord: slice(lat + window, lat - window),
                                    lon_coord: slice(lon - window, lon + window)})
    
        # remove pixels with area smaller or equal than the input pixel
        mask = upArea_window.where(upArea_window < area, np.nan)

        # from the remaining pixels, find that with the largest upstream area
        up_pixel = mask.where(mask == mask.max(), drop=True)

        return (up_pixel[lat_coord].round(4).item(),
                up_pixel[lon_coord].round(4).item())
                
    except Exception as e:
        print(f"Failed to find upstream pixel: {e}")
        return None, None

def downstream_pixel(
    lat: float,
    lon: float,
    upArea: xr.DataArray
) -> tuple:
    """
    This function finds the downstream pixel coordinates of a given point.
    
    Parameters:
    -----------
    lat: float
        Latitude of the input point.
    lon: float
        Longitude of the input point.
    upArea: xarray.DataArray
        Map of upstream area.
        
    Returns:
    --------
    tuple:
        Latitude and longitude of the immediate downstream pixel,
        or (None, None) if not found.
    """
    
    # Determine coordinate system
    lat_coord = 'lat' if 'lat' in upArea.coords else 'y'
    lon_coord = 'lon' if 'lon' in upArea.coords else 'x'
    
    try:
        # upstream area of the input coordinates
        area = upArea.sel({lat_coord: lat, lon_coord: lon}, method='nearest').values
        
        # spatial resolution of the input map
        resolution = np.mean(np.diff(upArea[lon_coord].values))
        
        # Define window around the input pixel
        window = 1.5 * resolution #np.array([-1.5 * resolution, 1.5 * resolution])
        upArea_window = upArea.sel({lat_coord: slice(lat + window, lat - window),
                                    lon_coord: slice(lon - window, lon + window)})
        
        # Mask out pixels with area less than or equal to the input pixel
        mask = upArea_window.where(upArea_window > area, drop=True)
        
        # Find the pixel with the smallest upstream area
        down_pixel = mask.where(mask == mask.max(), drop=True)
        
        return (down_pixel[lat_coord].round(4).item(),
                down_pixel[lon_coord].round(4).item())
                
    except Exception as e:
        print(f"Failed to find downstream pixel: {e}")
        return None, None

def outlet_width(
    chanbw: xr.DataArray,
    uparea: xr.DataArray,
    x: float,
    y: float,
    n_points: int = 5):
    """It searches for the width of the outlet of a water body. Given the coordinates of the lake/reservoir in the LISFLOOD grid, it finds the 'n_points' upstream and downstream, and it selects the width of the pixel with a largest change in the sequence.
    
    Parameters:
    -----------
    chanbw: xarray.DataArray
        static map of channel bankfull width
    uparea: xarray.DataArray
        static map of upstream area
    x: float
        coordinate X of the estimated outlet of the water body
    y: float
        coordinate Y of the estimated outlet of the water body
    n_points: int
        number of pixels to search. This number will be applied both upstream and downstream
        
    Returns:
    --------
    width: float
        width of the water body outlet
    """
    
    xs = np.array([x])
    ys = np.array([y])
    
    n_points = 2
    # coordinates of the two pixels upstream and downstream
    i = 0
    while i < n_points:
        try:
            y_up, x_up = upstream_pixel(ys[0], xs[0], uparea)
            xs = np.hstack([x_up, xs])
            ys = np.hstack([y_up, ys])
        except:
            break
        i += 1
        
    # coordinates of the two pixels downstream
    i = 0
    while i < n_points:
        try:
            y_down, x_down = downstream_pixel(ys[-1], xs[-1], uparea)
            if (y_down is not None) and (x_down is not None):
                xs = np.hstack([xs, x_down])
                ys = np.hstack([ys, y_down])
        except:
            break
        i += 1
        
    # extract channel width
    width = np.array([chanbw.sel(x=x, y=y, method='nearest').item() for x, y in zip(xs, ys)])
    
    return np.nanmin(width)


def duration_precip_indices(mask: pd.DataFrame) -> pd.Series:
    """
    Calculate the mean duration of precipitation extremes for each column in the mask DataFrame.

    This function assumes that the input DataFrame 'mask' contains boolean-like data (0s and 1s),
    where a 1 indicates the occurrence of a precipitation extreme and a 0 indicates no extreme.
    The function computes the mean duration of continuous sequences of 1s (extreme events) for
    each column, with the duration expressed in days.

    Parameters:
    -----------
    mask : pd.DataFrame
        A DataFrame with boolean-like values (0s and 1s), where the index represents time (as a
        DateTimeIndex) and each column corresponds to a different location or measurement.

    Returns:
    --------
    pd.Series
        A Series containing the mean duration of extreme precipitation events for each column in the mask
        DataFrame. The index of the Series corresponds to the columns of the input DataFrame.
    """
    
    duration = pd.Series(index=mask.columns, dtype=float)
    
    mask_diff = mask.astype(int).diff(axis=0)
    mask_diff.iloc[0] = mask.iloc[0].astype(int)
    for col in mask.columns:
        starts = mask_diff[mask_diff[col] == 1].index
        ends = mask_diff[mask_diff[col] == -1].index
        if len(starts) > len(ends):
            ends = ends.append(mask_diff.index[[-1]])
        elif len(starts) < len(ends):
            starts = mask_diff.index[[0]].append(starts)
        if len(starts) == 0:
            duration.loc[col] = 0
        else:
            duration.loc[col] = np.mean((ends - starts) / np.timedelta64(1, 'D'))
        
    return duration