import numpy as np
import pandas as pd
from typing import Union, List, Tuple, Dict, Optional
import logging
logger = logging.getLogger(__name__)

from .basemodel import Reservoir


class mHM(Reservoir):
    """Representation of reservoir routing in the Mesoscale Hydrological Model (mHM) as explained in Shrestha et al (2024)"""
    
    def __init__(
        self,
        Vmin: float,
        Vtot: float,
        Qmin: float,
        avg_inflow: float,
        avg_demand: float,
        w: float = 0.1, # Shin et al. (2019)
        alpha: float = 0.5, # called c*  in Shrestha et al. (2024). Default value from Hanasaki et al. (2006)
        beta: float = 1, # Shin et al. (2019)
        gamma: float = 0.85, # Shin et al. (2019)
        lambda_: float = 1,
        Atot: Optional[int] = None,
        timestep: int = 86400
    ):
        """        
        Parameters:
        -----------
        Vmin: float
            Volume (m3) associated to the conservative storage
        Vtot: float
            Total reservoir storage capacity (m3)
        Qmin: float
            Minimum outflow (m3/s)
        avg_inflow: float
            Average reservoir inflow (m3/s)
        avg_demand: float
            Average demand (m3/s)
        w: float
            Dimensionless parameter that controls the demand hedging
        alpha: float
            Dimensionless parameter that is a threshold that defines reservoirs whose releases are only based on demand (degree of regulation greater than alpha), or a combination of demand and inflow (otherwise)
        beta: float
            Dimensionless parameter that indirectly controls the proportion of inflow and demand in the releases
        gamma: float
            Dimensionless parameter that defines the normal storage: Vn = gamma * Vtot
        lambda_: float
            Dimensionless parameter that further controls the hedging in relation to the current reservoir filling
        Atot: integer (optional)
            Reservoir area (m2) at maximum capacity
        timestep: int
            Simulation time step in seconds.
        """
        
        assert 0 <= w <= 1, 'ERROR. Parameter "w" must be a value between 0 and 1'
        assert 0 <= alpha, 'ERROR. Parameter "alpha" (degree of regulation) must be positive'
        assert 0 <= gamma <= 1, 'ERROR. Parameter "gamma" must be a value between 0 and 1, as it represents the normal reservoir filling'
        
        # make sure that Vmin is not smaller than Vn
        Vmin = min(Vmin, gamma * Vtot)
            
        super().__init__(Vmin, Vtot, Qmin, Qf=None, Atot=Atot, timestep=timestep)
        
        # demand and degree of regulation
        self.avg_inflow = avg_inflow
        self.avg_demand = avg_demand
        self.dor = Vtot / (avg_inflow * 365 *24 * 3600) # called c in Shrestha et al. (2024)
        
        # reservoir parameters
        self.w = w
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.lambda_ = lambda_
        
        # normal storage
        self.Vn = gamma * Vtot
        # partition coefficient betweee demand-controlled (rho == 1) and non-demand-controlled reservoirs
        self.rho = min(1, (self.dor / alpha)**beta)
    
    def step(
        self,
        I: float,
        V: float,
        P: Optional[float] = None,
        E: Optional[float] = None,
        D: float = 0.0
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
        D: float
            Water demand (m3). It defaults to zero
            
        Returns:
        --------
        Q, V: List[float]
            Outflow (m3/s) and updated storage (m3)
        """
        
        # hedged demand
        exploitation = self.avg_demand / self.avg_inflow
        if exploitation >= 1 - self.w:
            hedged_demand = self.w * self.avg_inflow + (1 - self.w) * D / exploitation
        else:
            hedged_demand = D + self.avg_inflow - self.avg_demand             
            
        # outflow
        kappa = (V / self.Vn)**self.lambda_
        Q = self.rho * kappa * hedged_demand + (1 - self.rho) * I
        
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
            V = np.max([0, V - E * 1e-3 * A])
        # demand is not withdrawn here in this model, but later as a component of outflow
        
        # ouflow depending on the minimum outflow and storage level
        Q = np.max([self.Qmin, Q])
        eps = 1e-3
        if V - Q * self.timestep > self.Vtot:
            Q = (V - self.Vtot) / self.timestep + eps
        elif V - Q * self.timestep < self.Vmin:
            Q = (V - self.Vmin) / self.timestep - eps if V >= self.Vmin else 0
        if Q < 0:
            logger.warning(f'The simulated outflow was negative ({Q:.6f} m3/s). Limitted to 0')
            Q = 0
            
        # update reservoir storage with the inflow volume
        V -= Q * self.timestep
        
        assert 0 <= V, f'The volume at the end of the timestep is negative: {V:.0f} m3'
        assert V <= self.Vtot, f'The volume at the end of the timestep is larger than the total reservoir capacity: {V:.0f} m3 > {self.Vtot:.0f} m3'
            
        return Q, V
    
    def get_params(self):
        """It generates a dictionary with the reservoir paramenters in the model."""

        params = {
            'Vmin': self.Vmin,
            'Vn': self.Vn,
            'Vtot': self.Vtot,
            'Qmin': self.Qmin,
            'w': self.w,
            'alpha': self.alpha,
            'beta': self.beta,
            'gamma': self.gamma,
            'lambda': self.lambda_,
            'rho': self.rho,
            'Atot': self.Atot
        }
        # params = {key: float(value) for key, value in params.items()}
        params = {key: float(value) if value is not None else None for key, value in params.items()}

        return params