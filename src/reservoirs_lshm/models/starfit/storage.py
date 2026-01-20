import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from scipy.optimize import minimize
from pathlib import Path
from typing import Union, Optional, Dict, List, Literal
import logging
logger = logging.getLogger('storage')

from .inputs import read_reservoir_attributes, read_reservoir_data, rank_and_filter_data


def fit_storage(
    dam_id: int,
    storage_daily: Optional[pd.Series] = None,
    USRDATS_path: Optional[Union[str, Path]] = None,
    attributes: Optional[pd.DataFrame] = None,
    GRanD_path: Optional[Union[str, Path]] = None,
    capacity: str = 'CAP_MCM',
    cutoff_year: Optional[int] = None,
    min_days: int = 2920, # 8 years, the original repository uses 3650
    n_points: int = 3,
) -> Dict:
    """Fit parameters of storage targets
    
    Parameters:
    -----------
    dam_id: integer
        Dam ID in the GRanD database
    storage_daily: pandas.Series
        Daily time series of reservoir storage (hm3). The index must be the dates
    USRDATS_path: string or pathib.Path
        Path to the time series
    attributes: pandas.Series (optional)
        GRanD attributes for selected dam
    GRanD_path: string or pathlib.Path
        path to v1.3 of GRanD database. Only needed if 'attributes' is None
    capacity: string
        Field in the reservoir attributes used as reservoir storage capacity. By default "CAP_MCM"
    cutoff_year: integer (optional)
        Trim the time series to start this year
    min_days: integer
        Minimum number of days with storage values required to fit the target storage functions
    n_points: integer
        Number of maximum/minimum weekly storage values used to fit the flood/conservative storage harmonic function
    
    Returns:
    --------
    Dictionary
        id: integer
        weekly storage: pandas.DataFrame
            Weekly time series of median storage
        NOR upper bound: np.ndarray
            5 parameters of the flood storage harmonic function
        NOR lower bound: np.ndarray
            5 parameters of the conservative storage harmonic function
    """
    
    # extract reservoir storage capacity
    if attributes is None:
        attributes = read_reservoir_attributes(GRanD_path, dam_id)
    logger.info(f"Fitting targets for dam {dam_id}: {attributes['DAM_NAME']}")
    storage_capacity_MCM = attributes[capacity]
    
    if storage_daily is None:
        storage_daily = read_reservoir_data(USRDATS_path, dam_id).set_index('date')[['s_MCM']]
    else:
        if isinstance(storage_daily, pd.Series):
            storage_daily = pd.DataFrame(storage_daily)
            storage_daily.columns = ['s_MCM']
    storage_daily = storage_daily.loc[storage_daily.s_MCM.notnull()]
    
    if len(storage_daily) < min_days:
        logger.warning(f'Number of days with storage data ({len(storage_daily)}) smaller than the minimum ({min_days})')
        return {
            "id": dam_id,
            "capacity (MCM)": storage_capacity_MCM,
            "weekly storage": pd.DataFrame(),
            "NOR upper bound": [np.nan] * 5,
            "NOR lower bound": [np.nan] * 5
        }

    # # make sure that the timeseries has no missing days
    # start_date, end_date = storage_daily.date.min(), storage_daily.date.max()
    # storage_daily_filled = pd.DataFrame({'date': pd.date_range(start=start_date, end=end_date, freq='D')})
    # storage_daily = storage_daily_filled.merge(storage_daily, on='date', how='left')

    # check last year of data
    last_year_of_data = storage_daily.index.year.max()
    first_year_of_data = storage_daily.index.year.min()
    if cutoff_year is None or last_year_of_data < cutoff_year:
        cutoff_year = first_year_of_data
        logger.info(f"Cutoff year set back to {first_year_of_data}")

    # Convert to weekly storage (as % of capacity)
    storage_daily['year'] = storage_daily.index.year
    storage_daily['epiweek'] = storage_daily.index.isocalendar().week
    storage_daily = storage_daily[storage_daily.year >= cutoff_year]
    storage_weekly = (
        storage_daily.groupby(['year', 'epiweek'])
        # .agg(s_pct=('s_MCM', lambda x: round(100 * x.median() / storage_capacity_MCM, 2)))
        .agg(s_st=('s_MCM', lambda x: round(x.median() / storage_capacity_MCM, 2)))
        .reset_index()
    )
    storage_weekly = storage_weekly[storage_weekly['epiweek'].between(1, 52)]

    # Check for capacity and minimum violations
    # capacity_violations = storage_weekly[storage_weekly['s_pct'] > 100]
    capacity_violations = storage_weekly[storage_weekly['s_st'] > 1]
    minimum_violations = storage_weekly[storage_weekly['s_st'] < 0]
    if len(capacity_violations) > 0:
        logger.warning(f"{len(capacity_violations)} capacity violations found... ")
    if len(minimum_violations) > 0:
        logger.warning(f"{len(minimum_violations)} minimum violations found... ")

    # make sure that values don't exceed 0-100%
    # storage_weekly['s_pct'] = storage_weekly['s_pct'].clip(lower=0, upper=100)
    storage_weekly['s_st'] = storage_weekly['s_st'].clip(lower=0, upper=1)

    # For flood harmonic: rank the entries by descending 's_st' and keep the top 'n_points' for each epiweek
    data_for_flood_harmonic = rank_and_filter_data(storage_weekly, 's_st', n_points, ascending=False)
    
    # For conservation harmonic: rank the entries by ascending 's_st' and keep the top 'n_points' for each epiweek
    data_for_conservation_harmonic = rank_and_filter_data(storage_weekly, 's_st', n_points, ascending=True)

    # Fit the flood and conservation harmonics
    p_flood_harmonic = fit_constrained_harmonic(data_for_flood_harmonic).round(3)
    p_conservation_harmonic = fit_constrained_harmonic(data_for_conservation_harmonic).round(3)

    # Evaluate targets
    targets_flood = create_storage_harmonic(p_flood_harmonic, constrain=False)['target']
    targets_cons = create_storage_harmonic(p_conservation_harmonic, constrain=False)['target']

    # Adjust harmonic parameters based on targets
    if p_flood_harmonic[3] > targets_flood.max():
        p_flood_harmonic[3] = np.inf
    if p_flood_harmonic[4] < targets_flood.min():
        p_flood_harmonic[4] = -np.inf
    if p_conservation_harmonic[3] > targets_cons.max():
        p_conservation_harmonic[3] = np.inf
    if p_conservation_harmonic[4] < targets_cons.min():
        p_conservation_harmonic[4] = -np.inf

    return {
        "id": dam_id,
        "capacity (MCM)": storage_capacity_MCM,
        "weekly storage": storage_weekly,
        "NOR upper bound": p_flood_harmonic,
        "NOR lower bound": p_conservation_harmonic
    }


def fit_constrained_harmonic(data: pd.DataFrame) -> np.ndarray:
    """Fit parameters of a constrained harmonic function of target storage
    
    Parameters:
    -----------
    data: pandas.DataFrame
        Table with fields 'epiweek' and 's_st'
    
    Returns:
    --------
    pars: numpy.ndarray
        Fitted parameters of the constrained harmonic function of target storage
    """

    def evaluate_harmonic(pars: List):
        """Evaluate goodness-of-fit of fitted harmonic with the RMSE (root mean squared error). It is used as objective function in optimization of constrained harmonic.
        
        Parameters:
        -----------
        pars: list
            5 Parameters of the harmonic function of target storage
            
        Returns:
        --------
        rmse: float
            Root mean squared error
        """
        
        sin_term_vector = np.sin(2 * np.pi * data['epiweek'] / 52)
        cosin_term_vector = np.cos(2 * np.pi * data['epiweek'] / 52)
        fitted_harmonic = pars[0] + pars[1] * sin_term_vector + pars[2] * cosin_term_vector
        fitted_harmonic = np.minimum(np.maximum(fitted_harmonic, pars[4]), pars[3])
        
        rmse = np.sqrt(np.mean((data['s_st'] - fitted_harmonic)**2))
        return rmse
    
    # estimate the first 3 parameters by ordinary least squares
    initial_model = smf.ols('s_st ~ np.sin(2 * np.pi * epiweek / 52) + np.cos(2 * np.pi * epiweek / 52)', data=data).fit()
    intercept, sin_term, cosine_term = initial_model.params
    
    # estimate the last 2 parameters
    ub_on_curve = data['s_st'].quantile(0.9)
    lb_on_curve = data['s_st'].quantile(0.1)

    if (round(intercept, 5) == 100 or round(intercept, 5) == 0 or
       (round(sin_term, 5) == 0 and round(cosine_term, 5) == 0) or
       (round(ub_on_curve, 1) == round(lb_on_curve, 1))):
        return np.array([intercept, 0, 0, np.inf, -np.inf])

    optimized_constrained_harmonic = minimize(
        evaluate_harmonic,
        x0=[intercept, sin_term, cosine_term, ub_on_curve, lb_on_curve],
        bounds=[(0, None), (None, None), (None, None), (0, 100), (0, intercept)],
        method='L-BFGS-B'
    )
    pars = optimized_constrained_harmonic.x
    
    return pars


def create_storage_harmonic(
    parameters: List[float],
    freq: Literal['W', 'D'] = 'W',
    name: str = "target",
    constrain: bool = True
) -> pd.DataFrame:
    """It defines the weekly target values of storage given the parameters of the harmonic function
    
            storage = p1 + p2 · sin( 2 · pi · w · t ) + p3 · cos( 2 · pi · w · t )  # w: frequency, t: specific time
    
    Parameters:
    -----------
    parameters: numpy.ndarray
        vector of length 5 giving, in order, intercept, sine term, cosine term, and upper and lower constraints of the harmonic.
    freq: string
        Frequency of the harmonic function to be generated. Only two values are accepted: "W" for weekly, "D" for daily temporal resolution.
    name: string (optional)
        Character string naming the harmonic function. E.g., "flood" or "conservation." Default is simply "target"
    constrain: bool
        Constrain targets?
        
    Returns:
    --------
    storage_harmonic: pandas.Series
        The storage target levels by day/week of the year
    """
    
    # extract parameters
    p1, p2, p3 = parameters[:3]
    p4 = parameters[3] if constrain else float('inf')
    p5 = parameters[4] if constrain else float('-inf')

    # define harmonic function
    if freq == 'W':
        col, t = 'epiweek', 52
    elif freq == 'D':
        col, t = 'doy', 365
    else:
        raise ValueError(f'"freq" must be either "W" for weekly or "D" for daily, not "{freq}"')
    storage_harmonic = pd.DataFrame({col: np.arange(1, t + 1)})
    storage_harmonic[name] = (p1 +
                              p2 * np.sin(2 * np.pi * storage_harmonic[col] / t) +
                              p3 * np.cos(2 * np.pi * storage_harmonic[col] / t))
    storage_harmonic[name] = np.minimum(np.maximum(storage_harmonic[name], p5), p4)
    
    return storage_harmonic#.set_index('epiweek', drop=True)