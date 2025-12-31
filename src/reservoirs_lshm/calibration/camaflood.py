import numpy as np
import pandas as pd
from spotpy.objectivefunctions import kge
from spotpy.parameter import Uniform
from typing import List, Dict, Literal, Optional, Union
import logging
logger = logging.getLogger(__name__)

from .basecalibrator import Calibrator
from ..models import get_model
from .. import ParametersConfig
from ..utils.utils import return_period


class CamafloodCalibrator(Calibrator):
    """This class allows for calibrating 5 parameters in the Camaflood reservoir routine, 3 related to the storage limits, and 2 to the outflow limits.
    
    alpha: quantile of the storage records that defines the flood storage
            Vf = alpha * Vtot
    beta: defines the extreme storage as the distance between flood storage (Vf) and total capacity (Vtot)
            Ve = Vtot - beta * (Vtot - Vf)
    gamma: proportion of the flood storage (Vf) that corresponds to the normal storage (Vmin)
            Vmin = gamma * Vf
    delta: factor of the 100-year return period of inflow that defines the flood outflow (Qf)
            Qf = delta * Q100
    epsilon: factor of the mean inflow that defines the normal outflow (Qn)
            Qn = epsilon * Qf
    k: release coefficient that limits the outflow under flood conditions (V >= Vf & I > Qf)
    """
    
    def __init__(
        self,
        parameters: ParametersConfig,
        inflow: pd.Series,
        storage: pd.Series, 
        outflow: pd.Series, 
        Vmin: float, 
        Vtot: float, 
        catchment: int,
        Qmin: Optional[float] = None,
        precipitation: Optional[pd.Series] = None,
        evaporation: Optional[pd.Series] = None,
        demand: Optional[pd.Series] = None,
        Atot: Optional[float] = None,
        target: Union[Literal['storage', 'outflow'], List[Literal['storage', 'outflow']]] = 'storage', 
        obj_func=kge,
        spinup: Optional[int] = None
    ):
        """
        Parameters:
        -----------
        parametes: typed dictionary
            A dictionary that defines the parameters to be calibrated, together with the 'low' and 'high' values in
            the search range
        inflow: pd.Series
            Inflow time seris used to force the model
        storage: pd.Series
            Time series of reservoir storage
        outflow: pd.Series
            Observed outflow time series
        Vmin: float
            Volume (m3) associated to the conservative storage
        Vtot: float
            Total reservoir storage capacity (m3)
        catchment: integer
            Area (m2) of the reservoir catchment
        Qmin: float (optional)
            Minimum outflow (m3/s)
        precipitation: pandas.Series (optional)
            Time series of precipitation on the reservoir (mm)
        evaporation: pandas.Series (optional)
            Time series of open water evaporation from the reservoir (mm)
        demand: pandas.Series (optional)
            Time series of total water demand (m3)
        Atot: float (optional)
            Reservoir area (m2) at maximum capacity. Only needed if precipitaion or evaporation time series are provided as input
        target: list of strings
            Variable(s) targeted in the calibration. Possible values are 'storage' and/or 'outflow'
        obj_func:
            A function that assesses the performance of a simulation with a single float number. The optimization tries to minimize the objective function. We assume that the objective function would be either NSE or KGE, so the function is internally converted so that better performance corresponds to lower values of the objective function.
        spinup: integer (optional)
            Numer or time steps to use to warm up the model. These initial time steps will not be taken into account in the computation of model performance. By default, it is None and all the simulation will be used
        """
        
        super().__init__(inflow, storage, outflow, Vmin, Vtot, Qmin, precipitation, evaporation, demand, Atot, target, obj_func, spinup)
        
        self.catchment = catchment

        # define parameter distributions and names
        self.parameters = [
            Uniform(name=key, low=val['low'], high=val['high'])
            for key, val in parameters.items()
        ]
        self.param_names = list(parameters.keys())
        
    def pars2attrs(
        self,
        pars: List
    ) -> Dict:
        """It converts a list of model parameters into reservoir attributes to be used to declare a reservoir with `model.get_model()`
        
        Parameters:
        -----------
        pars: list
            Model parameters obtained, for instance, from the function `read_results()`

        Returns:
        --------
        attributes: dictionary
            Reservoir attributes needed to declare a reservoir using the function `models.get_model()`
        """

        # map parameter names and values
        param_dict = dict(zip(self.param_names, pars))
        
        # volume limits
        Vf = param_dict.get('alpha', 0.75) * self.Vtot 
        Ve = self.Vtot - param_dict.get('beta', 0.2) * (self.Vtot - Vf)
        Vmin = param_dict.get('gamma', 0.5) * Vf
        
        # outflow limits
        Qf = param_dict.get('delta', 0.3) * return_period(self.inflow, T=100)
        Qn = param_dict['epsilon'] * Qf if 'epsilon' in param_dict else self.inflow.mean()
        
        attributes = {
            'Vmin': Vmin, 
            'Vf': Vf,
            'Ve': Ve,
            'Vtot': self.Vtot,
            'Qn': Qn,
            'Qf': Qf,
            'catchment': self.catchment,
            'Atot': self.Atot
        }

        return attributes
        
    def simulation(
        self,
        pars: List[float],
    ) -> pd.Series:
        """Given a parameter set, it declares the reservoir and runs the simulation.
        
        Parameters:
        -----------
        pars: list of floats
            The set of parameter values to be simulated  
        
        Returns:
        --------
        sim: pandas.Series
            Simulated time series of the target variable
        """
        
        # declare the reservoir with the effect of the parameters
        reservoir_attrs = self.pars2attrs(pars)
        res = get_model('camaflood', **reservoir_attrs)
        if 'k' in self.param_names:
            res.k = dict(zip(self.param_names, pars))['k']
        self.reservoir = res
        
        # simulate
        sim = res.simulate(
            inflow=self.inflow, 
            Vo=self.observed['storage'].iloc[0],
            precipitation=self.precipitation,
            evaporation=self.evaporation,
            demand=self.demand
        )
        if self.spinup is not None:
            sim = sim.iloc[self.spinup:]
        
        return sim[self.target].round(2)