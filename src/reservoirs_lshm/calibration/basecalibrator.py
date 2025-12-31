import numpy as np
import pandas as pd
from spotpy.objectivefunctions import kge
from spotpy.parameter import Uniform
from typing import List, Literal, Optional, Dict, Union


class Calibrator(object):
    """Parent class used for the univariate calibration of reservoir modules. A specific child class needs to be created for each reservoir module to specify its parameter space and simulation process.
    """
    
    def __init__(
        self,
        inflow: pd.Series,
        storage: pd.Series, 
        outflow: pd.Series,
        Vmin: float, 
        Vtot: float, 
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
        inflow: pd.Series
            Inflow time seris used to force the model (m3/s)
        storage: pd.Series
            Observed storage time series (m3)
        outflow: pd.Series
            Observed outflow time series (m3/s)
        Vmin: float
            Volume (m3) associated to the conservative storage
        Vtot: float
            Total reservoir storage capacity (m3)
        Qmin: float
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
        
        # time series
        self.inflow = inflow
        self.precipitation = precipitation
        self.evaporation = evaporation
        self.demand = demand
        obs = pd.concat((outflow, storage), axis=1)
        obs.columns = ['outflow', 'storage']
        self.observed = obs
        
        # reservoir limits
        # volume
        self.Vmin, self.Vtot = Vmin, Vtot
        # outflow
        self.Qmin = Qmin
        # area
        self.Atot = Atot
        
        # target variable and objective function
        self.target = target
        self.obj_func = obj_func
        self.spinup = spinup
    
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

        pass
    
    def simulation(
        self,
        pars: List[float],
    ) -> pd.Series:
        """Given a parameter set, it declares the reservoir and runs the simulation.
        
        Parameters:
        -----------
        pars: List
            The set of parameter values to be simulated
            
        Returns:
        --------
        sim: pd.DataFrame
            Simulated time series of the target variable(s)
        """
        
        pass 

    def evaluation(
        self,
    ) -> pd.Series:
        """It extracts the observed time series of the target variable and removes (if necessary) the spinup time
            
        Returns:
        --------
        obs: pd.DataFrame
            Observed time series of the target variable
        """
        
        if self.spinup is not None:
            obs = self.observed[self.target].iloc[self.spinup:]
        else:
            obs = self.observed[self.target]
        
        return obs

    def objectivefunction(
        self,
        simulation: pd.Series,
        evaluation: pd.Series
    ) -> float:
        """It computes the objective function (self.obj_func) by comparison of the simulated and observed time series of the target variable 
        
        Parameters:
        -----------
        simulation: pd.Series
            Simulated time series
        evaluation: pd.Series
            Target time series
            
        Returns:
        --------
        of: float
            Value of the objective function
        """
        
        # compute the objective function
        of = []
        for var in self.target:
            perf = self.obj_func(evaluation[var], simulation[var])
            if isinstance(perf, tuple):
                of.append(1 - perf[0])
            else:
                of.append(1 - perf)
        
        return np.sqrt(np.sum(np.array(of)**2))