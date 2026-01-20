import numpy as np
import pandas as pd
from typing import Union, List, Tuple, Dict, Optional
import logging
logger = logging.getLogger(__name__)

from .basemodel import Reservoir


class Linear(Reservoir):
    """Representation of a linear reservoir"""
    
    def __init__(
        self,
        Vtot: float,
        T: int,
        Vmin: float = 0,
        Qmin: float = 0,
        Atot: Optional[int] = None,
        timestep: int = 86400
    ):
        """        
        Parameters:
        -----------
        Vtot: float
            Total reservoir storage capacity (m3)
        T: int
            Residence time in days. The coefficient of the linear reservoir is the inverse of T (1/T)
        Vmin: float
            Volume (m3) associated to the conservative storage
        Qmin: float
            Minimum outflow (m3/s)
        Atot: float (optional)
            Reservoir area (m2) at maximum capacity
        timestep: int
            Simulation time step in seconds.
        """
        
        super().__init__(Vmin=Vmin, Vtot=Vtot, Qmin=Qmin, Qf=None, Atot=Atot, timestep=timestep)
        
        # release coefficient
        self.k = 1 / (T * 86400) # s-1
        
    def step(
        self, 
        I: float, 
        V: float,
        P: Optional[float] = None,
        E: Optional[float] = None,
        D: Optional[float] = None,
    ) -> List[float]:
        """Given an inflow and an initial storage values, it computes the corresponding outflow
        
        Parameters:
        -----------
        I: float
            Inflow (m3/s)
        V: float
            Volume stored in the reservoir (m3)
        P: float (optional)
            Precipitaion on the reservoir (mm)
        E: float (optional)
            Open water evaporation (mm)
        D: float (optional)
            Consumptive demand (m3)
            
        Returns:
        --------
        Q, V: List[float]
            Outflow (m3/s) and updated storage (m3)
        """
        
        # estimate reservoir area at the beginning of the time step
        if P or E:
            if self.Atot:
                A = self.estimate_area(V)
            else:
                raise ValueError('To be able to model precipitation or evaporation, you must provide the maximum reservoir area ("Atot") in the reservoir declaration')
            
        # update reservoir storage
        V += I * self.timestep
        if P:
            V += P * 1e-3 * A
        if E:
            # evaporation can't happen if there's no water
            V = np.max([0., V - E * 1e-3 * A])
        if D:
            # demand can't withdraw water below the minimum storage
            V = np.max([self.Vmin, V - D])
        
        # outflow as a linear function of storage
        Q = np.max([self.Qmin, V * self.k])
        
        # limit outflow so the final storage is between 0 and 1
        # Q = np.max([np.min([Q, V / self.timestep]), (V - self.Vtot) / self.timestep])
        eps = 1e-3
        if V - Q * self.timestep > self.Vtot:
            Q = (V - self.Vtot) / self.timestep + eps
        elif V - Q * self.timestep < self.Vmin:
            Q = (V - self.Vmin) / self.timestep - eps if V >= self.Vmin else 0
        if Q < 0:
            logger.warning(f'The simulated outflow was negative ({Q:.6f} m3/s). Limitted to 0')
            Q = 0
            
        # update reservoir storage with the outflow volume
        V -= Q * self.timestep
        
        assert 0 <= V, f'The volume at the end of the timestep is negative: {V:.0f} m3'
        assert V <= self.Vtot, f'The volume at the end of the timestep is larger than the total reservoir capacity: {V:.0f} m3 > {self.Vtot:.0f} m3'
            
        return Q, V
    
    def get_params(self):
        """It generates a dictionary with the reservoir paramenters in the model."""

        params = {
            'Vmin': self.Vmin,
            'Vtot': self.Vtot,
            'Qmin': self.Qmin,
            'T': 1 / (self.k * 86400), # days
            'Atot': self.Atot
        }
        # params = {key: float(value) for key, value in params.items()}
        params = {key: float(value) if value is not None else None for key, value in params.items()}

        return params