from typing import Literal, Optional, Dict
import numpy as np
import pandas as pd
import logging
logger = logging.getLogger(__name__)

from .linear import Linear
from .lisflood import Lisflood
from .camaflood import Camaflood
from .mhm import mHM
from ..utils.utils import return_period

model_classes = {
    'linear': Linear,
    'lisflood': Lisflood,
    'camaflood': Camaflood,
    'mhm': mHM,
}

def get_model(
    model_name: Literal['linear', 'lisflood', 'camaflood', 'mhm'], 
    *args, 
    **kwargs
):
    """
    Creates an instance of the specified model class.
    
    Parameters:
    -----------
    model_name: string
        The name of the model class to instantiate.
    *args:
        Positional arguments to pass to the model class constructor.
    **kwargs:
        Keyword arguments to pass to the model class constructor.
        
    Returns:
    --------
    An instance of the specified model class.
    """
    
    # Get the class from the dictionary
    model_class = model_classes.get(model_name.lower())
    
    # Check if the model class exists
    if model_class is not None:
        # Create an instance of the model class
        return model_class(*args, **kwargs)
    else:
        raise ValueError(f"Model '{model_name}' not found. Available models: {list(model_classes.keys())}")
        
        
        
def default_attributes(
    model_name: Literal['linear', 'lisflood', 'camaflood', 'mhm'],
    inflow: pd.Series,
    Vtot: float,
    Vmin: Optional[float] = 0,
    Qmin: Optional[float] = 0,
    catchment: Optional[float] = None,
    Atot: Optional[float] = None,
    storage: Optional[pd.Series] = None,
    demand: Optional[pd.Series] = None,
) -> Dict:
    """It creates a dictionary with the model-specific attributes required to declare the reservoir class
    
    Parameters:
    -----------
    model_name: string
        Name of the reservoir model to be used: 'linear', 'lisflood', 'camaflood' or 'mhm'
    inflow: pandas.Series
        Time series of reservoir inflow (m3/s)
    Vtot: float (optional)
        Reservoir storage capacity (m3). Required by the 'linear', 'lisflood' and 'camaflood' models
    Vmin: float (optional)
        Minimum reservoir storage (m3). Required by the 'lisflood' model. If not provided, a value of 0 is used
    Qmin: float(optional)
        Minimum outflow (m3/s). Required by the 'lionear', 'lisflood' and 'camaflood' models. If not provided, a value of 0  is used
    catchment: float (optional)
        Reservoir catchment area (m2). Required by the 'camaflood' model
    storage: pandas.Series (optional)
        Time series of reservoir storage (m3). Required by the 'camaflood' routine
    demand: pandas.Series (optional)
        Time series of water demand (m3/s). Required by the 'mhm' routine
        
    Returns:
    --------
    attributes: dictionary
        Reservoir attributes needed to declare a reservoir using the function 'get_model'  in this same module
    """
          
    # model-generic dictionary of reservoir attributes
    attributes = {
        'Vtot': Vtot,
        'Vmin': Vmin,
        'Qmin': Qmin,
        'Atot': Atot if Atot is not None else None
    }
    
    # define model-specific reservoir attributes
    if model_name.lower() == 'linear':
        attributes.update({            
            'T': float(Vtot / (inflow.mean() * 24 * 3600))
        })
    elif model_name.lower() == 'lisflood':
        Vn, Vn_adj, Vf = np.array([0.67, 0.83, 0.97]) * Vtot
        Vmin = min(Vmin, Vn)
        Qf = .3 * return_period(inflow, T=100)
        Qn = min(inflow.mean(), Qf)
        attributes.update({
            'Vmin': Vmin,
            'Vn': Vn,
            'Vn_adj': Vn_adj,
            'Vf': Vf,
            'Qn': Qn,
            'Qf': Qf,
            'k': 1.2
        })
    elif model_name.lower() == 'camaflood':
        if catchment is None:
            raise ValueError(f"Model '{model_name}' requires the attribute 'catchment' reporting the reservoir catchment area (m2)")
        if storage is None:
            raise ValueError(f"Model '{model_name}' requires the attribute 'storage' reporting a time series of reservoir storage (m3)")
        Vf = float(storage.quantile(.75))
        Ve = Vtot - .2 * (Vtot - Vf)
        Vmin = .5 * Vf
        Qf = .3 * return_period(inflow, T=100)
        Qn = min(inflow.mean(), Qf)
        attributes.update({
            'Vf': Vf,
            'Ve': Ve,
            'Vmin': Vmin,
            'Qn': Qn,
            'Qf': Qf,
            'catchment': catchment
        })
        del attributes['Qmin']
    elif model_name.lower() == 'mhm':
        if demand is None:
            raise ValueError(f"Model '{model_name}' requires the attribute 'demand' reporting a time series of water demand (m3/s)")
        attributes.update({
            'gamma': 0.85, #float(storage.quantile(.9) / Vtot),
            'avg_inflow': inflow.mean(),
            'avg_demand': demand.mean()
        })
    else:
        raise ValueError('The model name must be one of the following: "linear", "lisflood", "camaflood", "mhm"')
        
    return attributes
    
    