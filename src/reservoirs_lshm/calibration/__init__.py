import pandas as pd
from typing import Literal, Optional, Union, Tuple, Dict, List
from pathlib import Path

from .linear import LinearCalibrator
from .lisflood import LisfloodCalibrator
from .camaflood import CamafloodCalibrator
from .mhm import MhmCalibrator


def get_calibrator(
    model_name: Literal['linear', 'lisflood', 'camaflood', 'mhm'], 
    *args, 
    **kwargs
):
    """
    Creates an instance of the specific calibration class for the reservoir model.
    
    Parameters:
    -----------
    model_name: string
        The name of the model class to instantiate. It must be one of the following values: 'linear', 'lisflood', 'camaflood' or 'mhm'
    *args:
        Positional arguments to pass to the calibrator class constructor.
    **kwargs:
        Keyword arguments to pass to the calibrator class constructor.
        
    Returns:
    --------
    An instance of the specified calibrator class.
    """
    
    if model_name.lower() == 'linear':
        return LinearCalibrator(*args, **kwargs)
    elif model_name.lower() == 'lisflood':
        return LisfloodCalibrator(*args, **kwargs)
    elif model_name.lower() == 'camaflood':
        return CamafloodCalibrator(*args, **kwargs)
    elif model_name.lower() == 'mhm':
        return MhmCalibrator(*args, **kwargs)
    else:
        raise ValueError("Invalid model name. Please choose either 'linear', 'lisflood', 'camaflood' or 'mhm'.")
        
        
def read_results(filename: Union[str, Path]) -> Tuple[pd.DataFrame, Dict]:
    """It reads the CSV file resulting from the calibration, and extracts the optimal parameter set.
    
    Parameters:
    -----------
    filename: string or pathlib.Path
        CSV file created during the calibration
        
    Returns:
    --------
    results: pandas.DataFrame
        The data in the CSV in Pandas format
    optimal_par: dictionary
        Optimal parameter set, i.e., those used in the iteration with lowest likelihood
    """
    
    # read results
    results = pd.read_csv(filename)
    results.index.name = 'iteration'
    parcols = {col: col[3:] for col in results.columns if col.startswith('par')}
    results.rename(columns=parcols, inplace=True)
    results['chain'] = results['chain'].astype(int)
    
    # extract optimal parameter set
    optimal_par = results.iloc[results.like1.idxmin()][parcols.values()].to_dict()
    
    return results, optimal_par