import os
os.environ['USE_PYGEOS'] = '0'
import pandas as pd
import geopandas as gpd
from typing import Union, Optional, List, Tuple, Dict
from pathlib import Path
import logging
logger = logging.getLogger('release')


def read_reservoir_data(
    USRDATS_path: Union[str, Path],
    dam_id: int
) -> pd.DataFrame:
    """Reads raw reservoir time series data for the specified dam

    Parameters:
    -----------
    USRDATS_path: string or pathlib.Path
        directory containing reservoir input time series
    dam_id: integer 
        id of dam; same as GRanD ID

    Returns:
    --------
    timeseries: pandas.DataFrame
        Daily time series of reservoir variables: 's_MCM' storage, 'i_cumecs' inflow, 'r_cumecs' release...
    """
    
    # read data
    USRDATS_path = Path(USRDATS_path) if isinstance(USRDATS_path, str) else USRDATS_path
    file_path = USRDATS_path / 'time_series_all' / f'ResOpsUS_{dam_id}.csv'
    if file_path.is_file():
        timeseries = pd.read_csv(
            file_path,
            usecols=['date', 'storage', 'inflow', 'outflow', 'elevation', 'evaporation'],
            parse_dates=['date']
        )
    else:
        logger.error(f'File not found: {file_path}')
        return None
    
    # make sure that the timeseries have no missing days
    start, end = timeseries.date.min(), timeseries.date.max()
    dates = pd.DataFrame({'date': pd.date_range(start=start, end=end, freq='D')})
    timeseries = dates.merge(timeseries, on='date', how='left')
    
    return timeseries.rename(columns={'storage': 's_MCM', 'inflow': 'i_cumecs', 'outflow': 'r_cumecs'})


def read_reservoir_attributes(
    GRanD_path: Union[str, Path],
    dam_id: Optional[int] = None
) -> pd.DataFrame:
    """Reads reservoir attributes from GRanD
    
    Parameters:
    -----------
    GRanD_path: string or pathlib.Path
        Directory containing the GRanD shapefile of reservoirs
    dam_id: integer (optional)
        Dam ID; same as GRanD ID. If None, attributes for all dams are returned

    Returns:
    --------
    attributes: pd.DataFrame
        Table of reservoir attributes for selected dams    
    """

    file_path = f"{GRanD_path}/GRanD_dams_v1_3.shp"
    if file_path.is_file():
        attributes = gpd.read_file(file_path)
        attributes.set_index('GRAND_ID', inplace=True, drop=False)
        # attributes_all = gdf[gdf['COUNTRY'] == "United States"].copy()
    else:
        logger.error(f'File not found: {file_path}')

    if dam_id is None:
        return attributes
    else:
        attributes = attributes[attributes.GRAND_ID == dam_id]
        assert len(attributes) == 1, "Dam ID should match exactly one dam."
        return attributes


# def read_GRanD_HUC8() -> pd.DataFrame:
#     """gets HUC8 for all US GRanD IDs
    
#     Returns:
#     --------
#     pandas.DataFrame of HUC8s
#     """
    
#     # Assuming that 'starfit' is the name of the directory where the 'extdata' folder is located
#     # and 'GRAND_HUC8.csv' is located inside the 'extdata' directory
#     file_path = "starfit/extdata/GRAND_HUC8.csv"
#     df = pd.read_csv(file_path, comment="#")
#     return df


def rank_and_filter_data(
    df: pd.DataFrame,
    rank_col: str,
    n_points: int,
    ascending: bool = True
) -> pd.DataFrame:
    """
    Rank the entries of the dataframe 'df' by the 'rank_col' and keep the top 'n_points'
    for each group defined by 'epiweek'. The ranking can be done in ascending or descending order.

    Parameters:
    -----------
    df: pandas DataFrame containing the data.
    rank_col: str
        the column on which to perform the ranking.
    n_points: int
        the number of top entries to keep for each group.
    ascending: bool, True for ascending order, False for descending order.

    Returns:
    --------
    A pandas DataFrame with the top 'n_points' ranked entries for each 'epiweek'.
    """
    ranked_data = (
        df.assign(
            rank=df.groupby('epiweek')[rank_col]
            .rank(ascending=ascending, method='first')
        )
        .query('rank <= @n_points')
        .drop(columns='rank')
        .sort_values('epiweek')
        .reset_index(drop=True)
    )
    ranked_data.epiweek = ranked_data.epiweek.astype(int)
    return ranked_data