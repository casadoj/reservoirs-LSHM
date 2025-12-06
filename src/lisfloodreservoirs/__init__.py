import yaml
from pathlib import Path
from typing import Union, Optional, List, Dict, TypedDict, Literal
import pandas as pd
from tqdm.auto import tqdm


class ParamRange(TypedDict):
    low: float
    high: float


ParameterName = Literal['alpha', 'beta', 'gamma', 'delta', 'epsilon', 'k']
ParametersConfig = Dict[ParameterName, ParamRange]


class Config:
    def __init__(self, config_file):
        
        # read configuration file
        with open(config_file, 'r', encoding='utf8') as ymlfile:
            self.cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)
        
        # data
        self.PATH_DATA = Path(self.cfg['data']['path'])
        self.RESERVOIRS_FILE = self.cfg['data']['reservoirs']
        self.PERIODS_FILE = self.cfg['data']['periods']
        path_results = Path(self.cfg['data'].get('results', './'))
        
        # model configuration
        self.MODEL = self.cfg['simulation']['model'].lower()
        self.PATH_DEF = path_results / f'{self.MODEL}' / 'default'
        self.PATH_DEF.mkdir(parents=True, exist_ok=True)
        self.INFLOW = self.cfg['simulation']['inputs'].get('inflow', 'inflow')
        self.PRECIPITATION = self.cfg['simulation']['inputs'].get('precipitation', None)
        self.EVAPORATION = self.cfg['simulation']['inputs'].get('evaporation', None)
        self.DEMAND = self.cfg['simulation']['inputs'].get('demand', None)
        self.SPINUP = self.cfg['simulation'].get('spinup', 0)
        
        # calibration
        cfg_cal = self.cfg.get('calibration', None)
        if cfg_cal is not None:
            # targets
            self.TARGET = cfg_cal.get('target', ['storage'])
            # parameters
            self.PARAMETERS: ParametersConfig = cfg_cal.get('parameters', None)
            # optimization
            self.MAX_ITER = cfg_cal['SCEUA'].get('max_iter', 2000)
            self.COMPLEXES = cfg_cal['SCEUA'].get('complexes', 8)
            self.KSTOP = cfg_cal['SCEUA'].get('kstop', 5)
            self.PEPS = cfg_cal['SCEUA'].get('peps', 0.01)
            self.PCENTO = cfg_cal['SCEUA'].get('pcento', 0.001)
            path_calib = path_results / self.MODEL / 'calibration'
            if len(self.TARGET) == 1:
                self.PATH_CALIB = path_calib / 'univariate' / self.TARGET[0]
            elif len(self.TARGET) == 2:
                self.PATH_CALIB = path_calib / 'bivariate'
            else:
                raise ValueError('ERROR. Only univariate or bivariate calibrations are supported')
            self.PATH_CALIB.mkdir(parents=True, exist_ok=True)
        
        
def read_attributes(
    path: Union[str, Path],
    reservoirs: Optional[List] = None,
    index_col: Optional[str] = 'GRAND_ID'
) -> pd.DataFrame:
    """It reads all the attribute tables from the specified dataset and, if provided, filters the selected reservoirs.
    
    Parameters:
    -----------
    path: string or pathlib.Path
        Directory where the dataset is stored
    reservoirs: list (optional)
        List of the reservoir ID selected
    index_col: string (otpional)
        Name of the column to be used as index
        
    Returns:
    --------
    attributes: pandas.DataFrame
        Concatenation of all the attributes in the dataset
    """      
        
    # import all tables of attributes
    try:
        attributes = pd.concat([pd.read_csv(file, index_col=index_col) for file in path.glob('*.csv')],
                               axis=1,
                               join='outer')
        attributes.index.name = index_col
        if reservoirs is not None:
            if isinstance(reservoirs, list) is False:
                reservoirs = [reservoirs]
            attributes = attributes.loc[reservoirs]
    except Exception as e:
        raise ValueError(f'ERROR while reading attribute tables from directory {path}: {e}') from e
        
    return attributes


def read_timeseries(
    path: Union[str, Path],
    reservoirs: Optional[List[int]] = None,
    periods: Optional[Dict[int, Dict[str, pd.Timestamp]]] = None,
    variables: Optional[List[str]] = None,
) -> Dict[int, pd.DataFrame]:
    """It reads the time series in the dataset and saves them in a dictionary.
    
    Parameters:
    -----------
    path: string or pathlib.Path
        Directory where the dataset is stored
    reservoirs: list (optional)
        List of the reservoir ID selected
    periods: dictionary (optional)
        If provided, it cuts the time series to the specified period. It is a dictionary of dictionaries, where the keys are the reservoir ID, and the values are dictionaries with two entries ('start' and 'end') that contain timestamps of the selected beginning and end of the study period
    variables: list (optional)
        If selects the time series to be loaded. If not provided, it searches for columns 'inflow', 'storage', 'outflow', 'elevation'
    
    Returns:
    --------
    timeseries: dictionary
        It contains the timeseries of the selected reservoirs as pandas.DataFrame
    """

    if reservoirs is None:
        reservoirs = [int(file.stem) for file in path.glob('*.csv')]

    # if variables is None:
    #     variables = ['inflow', 'storage', 'outflow', 'elevation']
        
    # read time series
    timeseries = {}
    for ID in tqdm(reservoirs):
        # read time series
        file = path / f'{ID}.csv'
        if file.is_file():
            ts = pd.read_csv(file)
            ts['date'] = pd.to_datetime(ts['date'])
            ts = ts.set_index('date').sort_index().asfreq('D')
        else:
            print(f"File {file} doesn't exist")
            continue

        # select study period
        try:
            if periods is not None:
                start, end = [periods[str(ID)][f'{x}_dates'][0] for x in ['start', 'end']]
                ts = ts.loc[start:end, :]
        except Exception as e:
            print(f'Error while trimming to the study period the time series for ID {ID}:\m{e}')

        # select varibles
        try:
            if variables is not None:
                missing_vars = set(variables).difference(ts.columns)
                if len(missing_vars) > 0:
                    print(f'Time series for ID {ID} is missing variables: {missing_vars}')
                ts = ts[ts.columns.intersection(variables)]
            # convert storage variables to m3
            ts.iloc[:, ts.columns.str.contains('storage')] *= 1e6
        except Exception as e:
            print(f'Error while selecting variables from the time series for ID {ID}:\m{e}')

        # save time series
        try:
            timeseries[ID] = ts
        except Exception as e:
            print(f'Time series for ID {ID} could not be saved:\n{e}')
        
    return timeseries