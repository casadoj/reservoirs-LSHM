import numpy as np
import pandas as pd
import xarray as xr
from typing import Union, List, Dict, Tuple, Optional, Literal
from statsmodels.distributions.empirical_distribution import ECDF
from timezonefinder import TimezoneFinder


def quantile_mapping(obs: pd.Series, sim: pd.Series) -> pd.Series:
    """It corrects the bias in the "sim" time series to replicate the empirical cumulative density function of the observed time series
    
    Parameters:
    -----------
    obs: pandas.Series
        Observed time series. May contain missing values
    sim: pandas.Series
        Simulated time series. May not contain missing values
        
    Returns:
    --------
    sim_bc: pandas.Series
        Bias corrected time series of equal length as "sim"
    """
    
    # remove timesteps with missing values
    data = pd.concat((obs, sim), axis=1)
    data.columns = ['obs', 'sim']
    data.dropna(axis=0, how='any', inplace=True)

    # empirical cumulative density functions
    obs_ecdf = ECDF(data.obs)
    sim_ecdf = ECDF(data.sim)

    # quantiles of the input "sim"
    sim_quantiles = sim_ecdf(sim)

    # bias correct "sim"
    x = data.obs.sort_values().unique()
    sim_bc = np.interp(sim_quantiles, obs_ecdf(x), x)
    
    return pd.Series(sim_bc, index=sim.index)


# def decomposition(data: Union[pd.DataFrame, pd.Series]) -> Tuple[Union[pd.DataFrame, pd.Series], Union[pd.DataFrame, pd.Series], Union[pd.DataFrame, pd.Series]]:
#     """It decomposes the timeseries in three components: annual average, seasonality and residuals.
    
#     Parameters:
#     -----------
#     data: Union[pd.DataFrame, pd.Series]
#         Time series to be decomposed
        
#     Returns:
#     --------
#     annual: Union[pd.DataFrame, pd.Series]
#         Timeseries of mean annual values. The length of the timeseries will be the number of years in "data"
#     seasonality: Union[pd.DataFrame, pd.Series]
#         Timeseries of montly mean values after removal of the annual trend. The length of the timeseries will be alwas 12, the months
#     residuals: Union[pd.DataFrame, pd.Series]
#         Timeseries of residuals, that is, the difference of the original "data" and the annual and seosonal timeseris. The length of the timeseries is the same as the input "data"
#     """
    
#     # annual average storage
#     annual = data.resample('Y').mean()
#     annual.index = annual.index.year

#     # seasonal variability
#     detrended = data - data.resample('Y').transform('mean')
#     seasonality = detrended.groupby(detrended.index.month).mean()

#     # residual
#     residuals = detrended - detrended.groupby(detrended.index.month).transform('mean')
    
#     return annual, seasonality, residuals


# class Decomposition:
#     def __init__(self, original: Tuple[Union[pd.DataFrame, pd.Series]], trend: Tuple[Union[pd.DataFrame, pd.Series]], seasonal: Tuple[Union[pd.DataFrame, pd.Series]], residuals: Tuple[Union[pd.DataFrame, pd.Series]]):
#         self.original = original
#         self.trend = trend
#         self.seasonal = seasonal
#         self.residual = residuals

        
        
# def decompose_timeseries(data: Union[pd.DataFrame, pd.Series], window: int = 365, center: bool = True) -> Decomposition:
#     """It decomposes the timeseries in three components: trend, seasonality and residuals.
    
#     Parameters:
#     -----------
#     data: Union[pd.DataFrame, pd.Series]
#         Time series to be decomposed
        
#     Returns:
#     --------
#     DecompositionResult:
#         Object with three methods: trend(), seasonal(), residuals()
#     """ 

#     # trend as the 365 rolling mean
#     trend = data.rolling(window=365, min_periods=180, center=center).mean()

#     # seasonality
#     detrended = data - trend
#     seasonal = detrended.groupby(detrended.index.month).transform('mean')

#     # residuals
#     residual = detrended - seasonal

#     return Decomposition(data, trend, seasonal, residual)


def decompose_timeseries(
    data: pd.Series, 
    window: int = 365, 
    center: bool = True
) -> pd.DataFrame:
    """It decomposes the timeseries in three components: trend, seasonality and residuals.
    
    Parameters:
    -----------
    data: pandas.Series
        Time series to be decomposed
    window: integer
        Width of the rolling window used to extract the trend
    center: boolean
        Whether to center or not the rolling window used to extract the trend        
        
    Returns:
    --------
    DecompositionResult:
        Object with three methods: trend(), seasonal(), residuals()
    """ 
    
    assert isinstance(data, pd.Series), '"data" must be a pandas.Series'
    
    # trend as the 365 rolling mean
    trend = data.rolling(window=window, min_periods=180, center=center).mean()

    # seasonality
    detrended = data - trend
    seasonal = detrended.groupby(detrended.index.month).transform('mean')

    # residuals
    residual = detrended - seasonal
    
    decomposition = pd.concat((data, trend, seasonal, residual), axis=1)
    decomposition.columns = ['original', 'trend', 'seasonal', 'residual']
    
    return decomposition


def clean_storage(
    storage: pd.Series,
    w: int = 7,
    error_thr: float = .1,
    inplace: bool = False
    ) -> Optional[pd.Series]:
    """It removes values from a reservoir storage time series that exceed an acceptable error ("error_thr") between the actual value and that of a centered, moving median of width "w"
    
    Parameters:
    -----------
    storage: pandas.Series
        Time series of reservoir storage
    w: integer
        Width of the window used for the moving median
    error_thr: float
        Maximum acceptable error. Time steps with an absolute error above this threshold will be converted into NaN
    inplace: boolean
        Whether to affect the input time series (True), or create a new corrected time series (False)
    
    Returns:
    --------
    storage_: pandas.Series (optional)
        If "inplace" is False, the function returns the corrected time series
    """
    
    storage_ = storage.copy()
    
    # remove values lower or equal than 0
    storage_[storage <= 0] = np.nan
    
    median = storage_.rolling(w, center=True, min_periods=int(np.floor(w / 2))).median()
    error = (storage_ - median) / median
    storage_[error.abs() > error_thr] = np.nan
    
    if inplace:
        storage = storage_
    else:
        return storage_

    
def clean_storage2(
    storage: pd.Series,
    capacity: Optional[float] = None,
    w: int = 7,
    error_thr: float = .1,
    inplace: bool = False
    ) -> Optional[pd.Series]:
    """It removes values from a reservoir storage time series that exceed an acceptable error ("error_thr") between the actual value and that of a centered, moving median of width "w". Compared with the previous version of this function, the error is a quotient of the storage capacity, instead of the median storage
    
    Parameters:
    -----------
    storage: pandas.Series
        Time series of reservoir storage
    capacity: float (optional)
        Storage capacity of the reservoir. If not provided, it is assumed as the maximum of the "storage" time series (susceptible to outliers)
    w: integer
        Width of the window used for the moving median
    error_thr: float
        Maximum acceptable error. Time steps with an absolute error above this threshold will be converted into NaN
    inplace: boolean
        Whether to affect the input time series (True), or create a new corrected time series (False)
    
    Returns:
    --------
    storage_: pandas.Series (optional)
        If "inplace" is False, the function returns the corrected time series
    """
    
    storage_ = storage.copy()

    if capacity is None:
        capacity = storage.max()
            
    # remove values lower or equal than 0
    storage_[storage == 0] = np.nan
        
    # remove outliers
    median = storage_.rolling(w, center=True, min_periods=int(np.floor(w / 2))).median()
    error = (storage_ - median) / capacity # only change compared with `clean_storage2`
    storage_[error.abs() > error_thr] = np.nan
    
    if inplace:
        storage = storage_
    else:
        return storage_
    
    
def clean_inflow(
    inflow: pd.Series,
    storage: Optional[pd.Series] = None,
    outflow: Optional[pd.Series] = None,
    grad_thr: int = 10000,
    balance_thr: Optional[int] = 5,
    int_method: Literal['linear', 'quadratic', 'spline'] = 'linear',
    inplace: bool = False,
    **kwargs
    ) -> pd.Series:
    """It prepares the time series of inflow to be used as input for the reservoir model, i.e., it removes erroneous peaks and it fills in the gaps.
    
        - To remove peaks, it uses two error identifiers. When both conditions identify an error in the inflow time series, that time step is converted into a NaN:
            1. If a storage and an outflow time series are provided, we can estimate the mass balance error. It this error exceeds an acceptable value ("balance_thr"), that time step is identified as a possible inflow error.
            2. The gradient (the difference between the current and the previous inflow). Time steps were the absolute gradient exceeds a maximum acceptable value ("grad_thr") are identified as possible inflow errors.
    
        - Three Pandas standard methods are available to fill in the gaps ('linear', 'quadratic', 'spline').
    
    Parameters:
    -----------
    inflow: pandas.Series
        Time series of reservoir inflow
    storage: pandas.Series (optional)
        Time series of reservoir storage. This time series must have been clean before, as it will be use as reference in the mass balance analysis
    outflow: pandas.Series (optional)
        Time series of reservoir outflow.
    grad_thr: integer
        Maximum acceptable streamflow gradient. The units are those of the inflow time series, usually m3/s. Time steps with a gradient above this threshold will be converted into NaN
    balance_thr: integer
        Maximum acceptable error in the mass balance. Dimensionless, it represents the relative error between mass-balance estimated storage and reported storage. It is only needed if "storage" and "outflow" are provided
    int_method: string
        Interpolation method used to fill in gaps in the time series
    inplace: boolean
        Whether to overwrite the input time series (True), or create a new corrected time series (False)
        
    Keyword arguments:
    ------------------
    limit: integer
        Maximum number of consecutive NaNs to fill. Must be greater than 0       
    
    Returns:
    --------
    storage_: pandas.Series (optional)
        If "inplace" is False, the function returns the corrected time series
    """
    
    inflow_ = inflow.copy()
    
    # remove values higher than 100,000 m3/s
    inflow_[inflow >= 1e5] = np.nan
    
    if (storage is not None) and (outflow is not None):
        # temporal resolution in seconds
        At = inflow.index.to_series().diff().mode()[0].total_seconds()
        # estimate storage based on the mass balance: S1 = S0 + (I1 - O1) * At
        estimated_storage = storage.shift(1) + (inflow - outflow) * At * 1e-6
        # error between reported and estimated storage
        error = (storage - estimated_storage) / storage
        clean_inflow.error = error
        # time steps that exceed the acceptable error
        mask_error = error.abs() > balance_thr
    else:
        mask_error = True
    
    # time steps that exceed the acceptable gradient
    grad = inflow_.diff()
    clean_inflow.gradient = grad
    mask_grad = grad.abs() > grad_thr
    
    # remove values that exceed the threshold(s)
    inflow_[mask_error & mask_grad] = np.nan
    
    # fill in gaps
    inflow_ = inflow_.interpolate(method=int_method, limit=kwargs.get('limit', 7), order=2)
    
    # make sure that the filling didn't produce negative values ('quadratic' or 'spline')
    inflow_[inflow_ < 0] = 0
    
    if inplace:
        inflow = inflow_
    else:
        return inflow_
    
    
def create_demand(
    outflow: pd.Series,
    water_stress: float = 1,
    window: int = 8
    ) -> pd.Series:
    """It creates a demand time series out of the outflow time series.
    First, it computes the average outflow for every day of the year (scaled with the "water_stress" coefficient), and then it applies a rolling mean of "window" width to smooth the resulting time series.
    
    Parameters:
    -----------
    outflow: pandas.Series
        Observed outflow time series
    water_stress: float
        Coefficient between average demand and average inflow
    window: integer
        Width of the rolling window used to smooth the demand time series
    
    Returns:
    --------
    demand: pandas.Series
        Time series of water demand
    """
    
    assert 0 < water_stress <= 1, "'water_stress' must be a value in the range (0, 1]"
    assert window > 0, "'window' must be an integer larger than 0"
    
    # average demand of every day of the year
    demand_d = outflow.groupby([outflow.index.month, outflow.index.day]).mean() * water_stress
    
    # smooth demand by applying a rolling mean
    buffer = int(window / 2)
    demand_d = pd.concat((demand_d.iloc[-buffer:], demand_d, demand_d.iloc[:buffer - 1]), axis=0)
    demand_d = demand_d.rolling(window, center=True).mean()
    demand_d.dropna(inplace=True)
    demand_d.index = np.arange(1, demand_d.shape[0] + 1)
    
    # define demand time series for the complete simulation period
    demand = pd.Series({date: demand_d.loc[date.dayofyear] for date in outflow.index})
    
    return demand


def define_period(series: Union[pd.Series, pd.DataFrame]) -> Tuple[np.datetime64, np.datetime64]:
    """
    If finds the beginning and end of the longest period with available data. If that period does not exceed a minimum number of years, the function skips the computation
    
    Parameters:
    -----------
    series: pandas.Series or pandas.DataFrame
        Time series to be checked 
    
    Returns:
    --------
    start: numpy.datetime64
        Beginning of the longest period
    end: numpy.datetime64
        End of the longest period
    """

    # boolean series of available data
    available = ~series.isnull()
    if isinstance(available, pd.DataFrame):
        available = available.all(axis=1)
    
    # find starts and ends of periods of years with enough inflow data
    diff = available.astype(int).diff()
    diff.iloc[0] = available.iloc[0].astype(int)
    starts = diff[diff == 1].index.values
    ends = diff[diff == -1].index.values
    if len(starts) == len(ends) + 1:
        ends = np.hstack((ends, diff.index.values[-1]))

    # return start and end of the longest period
    durations = (ends - starts) / np.timedelta64(1, 'D')
    if len(durations) > 0:
        dmax, imax = durations.max(), durations.argmax()
        return starts[imax], ends[imax] - np.timedelta64(1, 'D')
    else:
        return np.nan, np.nan
    
    
def time_encoding(
    da: xr.DataArray,
    period: int
) -> Tuple[xr.DataArray, xr.DataArray]:
    """
    Transforms time feature values in an xarray.DataArray to sine and cosine components.

    Parameters:
    -----------
    da: xarray.DataArray)
        An xarray.DataArray with time feature values (e.g., month, day of year).
    period: integer
        The period of the time feature (e.g., 12 for months, 7 for days of the week).

    Returns:
    --------
    sin_da, cos_da (tuple of xarray.DataArray):
        Sine and cosine transformations of the time feature values.
    """
    
    # Normalize time feature values to [0, 2π]
    if da.min().item() == 1:
        norm_da = (da - 1) * 2 * np.pi / period
    elif da.min().item() == 0:
        norm_da = da * 2 * np.pi / period
    else:
        norm_da = (da - 1) * 2 * np.pi / period
        
    # correct leap years, if necessary
    norm_da = norm_da.where(norm_da <= np.pi * 2, np.pi * 2)
    
    return np.sin(norm_da), np.cos(norm_da)