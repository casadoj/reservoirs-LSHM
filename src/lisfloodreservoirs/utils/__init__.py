import yaml
from pathlib import Path
from typing import Union, Optional, List, Dict
# from datetime import datetime
import pandas as pd

class DatasetConfig:
    
    def __init__(self, config_file: Union[str, Path]):
        
        # read configuration file
        with open(config_file, 'r', encoding='utf8') as ymlfile:
            self.cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)
        
        self.VERSION = self.cfg['version']
        
        # input paths
        self.PATH_LISFLOOD = Path(self.cfg['paths']['LISFLOOD']['root'])
        self.PATH_RESOPS = Path(self.cfg['paths']['ResOps']['root'])
        self.PATH_OBS_TS =  self.PATH_RESOPS / self.cfg['paths']['ResOps']['obs_timeseries']
        self.PATH_SIM_TS = self.PATH_RESOPS / self.cfg['paths']['ResOps']['sim_timeseries']
        self.PATH_GRAND = Path(self.cfg['paths']['GRanD'])
        
        # output paths
        self.PATH_ATTRS = self.PATH_RESOPS / self.VERSION / 'attributes'
        self.PATH_ATTRS.mkdir(parents=True, exist_ok=True)
        self.PATH_TS = self.PATH_RESOPS / self.VERSION / 'time_series'
        self.PATH_TS.mkdir(parents=True, exist_ok=True)
        
        # period
        self.START = pd.to_datetime(self.cfg['period'].get('start', '1900-01-01'))
        self.END = pd.to_datetime(self.cfg['period'].get('end', 'now')) #.date()
        
        self.NORMALIZE = self.cfg['normalize']
        
        # conditions
        self.MIN_AREA = self.cfg['conditions']['min_area']
        self.MIN_VOL = self.cfg['conditions']['min_volume']
        self.MIN_DOR = self.cfg['conditions']['min_dor']
        self.MIN_DOD = self.cfg['conditions']['min_dod']
        self.MIN_YEARS = self.cfg['conditions']['min_years']
        self.TOL_BIAS = self.cfg['conditions']['tol_bias']
        

class APIConfig:
    
    def __init__(self, config_file: Union[str, Path]):
        
        with open(config_file, 'r', encoding='utf8') as ymlfile:
            self.cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)
            
        self.URL = self.cfg['url']
        self.USERNAME = self.cfg['username']
        self.PASSWORD = self.cfg['password']