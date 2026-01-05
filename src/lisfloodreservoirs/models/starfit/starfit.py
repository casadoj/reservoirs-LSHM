import numpy as np
import pandas as pd
from typing import List, Optional, Literal
from tqdm.auto import tqdm
import logging
logger = logging.getLogger(__name__)

from .storage import create_storage_harmonic
from .release import create_release_harmonic
from ..basemodel import Reservoir


class Starfit(Reservoir):
    """
    Starfit is a subclass of the Reservoir class that models a reservoir with specific
    harmonic storage and release patterns for flood and conservation purposes.

    Parameters:
    -----------
    Vtot: float
        The total volume of the reservoir [MCM].
    avg_inflow (float): 
        The average inflow into the reservoir [MCM/day].
    pars_Vf: List 
        Parameters defining the harmonic storage pattern for flood conditions.
    pars_Vc: List 
        Parameters defining the harmonic storage pattern for conservation.
    pars_Qharm: List 
        Parameters defining the harmonic release pattern from the reservoir.
    pars_Qresid: List 
        Parameters for calculating residual releases from the reservoir.
    Qmin: float 
        The minimum allowable release from the reservoir [MCM/day].
    Qmax: float 
        The maximum allowable release from the reservoir [MCM/day].
    Atot: integer (optional)
        Reservoir area (m2) at maximum capacity

    Attributes:
    -----------
    avg_inflow: float 
        Stores the average inflow value provided during initialization.
    NOR: pandas.DataFrame
        A pandas DataFrame containing the normalized operational rules for flood and conservation storage, indexed by day of the year.
    Qharm: pandas.Series
        A pandas Series containing the harmonic release pattern, indexed by day of the year.
    parsQresid: List
        Stores the parameters for residual releases.
    Qmax: float 
        The maximum allowable release from the reservoir.

    Methods:
    --------
    Inherits all methods from the Reservoir class and does not define any new explicit methods.

    Notes:
    ------
    The class extends the functionality of the Reservoir base class by incorporating
    additional attributes related to harmonic storage and release patterns. It uses
    daily frequencies for these patterns and sets up the operational rules based on
    the input parameters.
    """
    
    def __init__(
        self,
        Vtot: float,
        avg_inflow: float,
        pars_Vf: List,
        pars_Vc: List,
        pars_Qharm: List,
        pars_Qresid: List,
        Qmin: float,
        Qmax: float,
        Atot: Optional[int] = None,
    ):
        
        super().__init__(Vmin=None, Vtot=Vtot, Qmin=Qmin, Qf=None, Atot=Atot, timestep=86400)
        
        # self.Vtot = Vtot
        self.avg_inflow = avg_inflow
        self.NOR = pd.concat((create_storage_harmonic(pars_Vf, freq='D', name='flood').set_index('doy'),
                              create_storage_harmonic(pars_Vc, freq='D', name='conservation').set_index('doy')),
                             axis=1)
        self.Qharm = create_release_harmonic(pars_Qharm, freq='D').set_index('doy').squeeze()
        self.parsQresid = pars_Qresid
        # self.Qmin = Qmin
        self.Qmax = Qmax
        
    def step(
        self, 
        I: float,
        V: float,
        doy: int,
        P: Optional[float] = None,
        E: Optional[float] = None,
        D: Optional[float] = None,
    ) -> List[float]:
        """Given an inflow and an initial storage values, it computes the corresponding outflow and storage at the end of the timestep
        
        Parameters:
        -----------
        I: float
            Inflow (m3/s)
        V: float
            Volume stored in the reservoir (m3)
        doy: integer
            Doy of the year. It must be a value between 1 and 365
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
            V = np.max([0, V - E * 1e-3 * A])
        if D:
            # demand can't withdraw water below the minimum storage
            # V = max(self.Vmin, V - D) # Starfit doesn't define Vmin
            V = np.max([0, V - D])
        
        # standardised inputs
        I_st = I / self.avg_inflow - 1
        V_st = V / self.Vtot
        
        # flood/conservation storage that week
        doy = 365 if doy == 366 else doy
        assert 1 <= doy <= 365, f'"doy" must be a value between 1 and 365 (including both): {doy} was provided'
        Vf, Vc = self.NOR.loc[doy, ['flood', 'conservation']]            
        
        # harmonic component of the release
        harm = self.Qharm[doy]
        # residual component of the release
        A_t = (V_st - Vc) / (Vf - Vc) # storage availability
        eps = self.parsQresid['Intercept'] + A_t * self.parsQresid['a_st'] + I_st * self.parsQresid['i_st']      
        # normal release
        Qnor = self.avg_inflow * (harm + eps + 1)
        
        # compute release
        if V_st < Vc:
            # Q = self.Qmin # original routine
            Q = min(self.Qmin + (Qnor - self.Qmin) * V_st / Vc, I)
        elif Vc <= V_st <= Vf:
            # Q = min(Qnor, self.Qmax) # original routine
            Q = max(min(Qnor, self.Qmax), self.Qmin)
        elif V_st > Vf:
            # Q = min((V_st - Vf) * self.Vtot / self.timestep + I, self.Qmax) # original routine
            Q = Qnor + (self.Qmax - Qnor) * (V_st - Vf) / (1 - Vf)

        # ensure mass conservation
        # Q = max(min(Q, V / self.timestep), (V - self.Vtot) / self.timestep)
        eps = 1e-3
        if V - Q * self.timestep > self.Vtot:
            Q = (V - self.Vtot) / self.timestep + eps
        elif V - Q * self.timestep < 0:
            Q = V / self.timestep - eps
        if Q < 0:
            logger.warning(f'The simulated outflow was negative ({Q:.6f} m3/s). Limitted to 0')
            Q = 0
        
        # update storage
        V -= Q * self.timestep
        
        assert 0 <= V, f'The volume at the end of the timestep is negative: {V:.0f} m3'
        assert V <= self.Vtot, f'The volume at the end of the timestep is larger than the total reservoir capacity: {V:.0f} m3 > {self.Vtot:.0f} m3'
        assert 0 <= Q, f'The simulated outflow is negative: {Q:.6f} m3/s'
            
        return Q, V
    
    def simulate(
        self,
        inflow: pd.Series,
        Vo: Optional[float ] = None,
        precipitation: Optional[pd.Series] = None,
        evaporation: Optional[pd.Series] = None,
        demand: Optional[pd.Series] = None,
    ) -> pd.DataFrame:
        """Given an inflow time series (m3/s) and an initial storage (m3), it computes the time series of outflow (m3/s) and storage (m3)
        
        Parameters:
        -----------
        inflow: pd.Series
            Time series of flow coming into the reservoir (m3/s)
        Vo: float (optional)
            Initial value of reservoir storage (m3). If not provided, it is assumed that the normal storage is the initial condition
        demand: pandas.Series (optional)
            Time series of total water demand
            
        Returns:
        --------
        pd.DataFrame
            A table that concatenates the storage (m3), inflow (m3/s) and outflow (m3/s) time series.
        """
        
        if Vo is None:
            Vo = .5 * self.Vtot
        
        if precipitation is not None and not isinstance(precipitation, pd.Series):
            raise ValueError('"precipitation" must be a pandas.Series representing a time series of precipitation (mm) on the reservoir.')
        if evaporation is not None and not isinstance(evaporation, pd.Series):
            raise ValueError('"evaporation" must be a pandas.Series representing a time series of open water evaporation (mm) from the reservoir.')
        if demand is not None and not isinstance(demand, pd.Series):
            raise ValueError('"demand" must be a pandas.Series representing a time series of water demand (m3/s).')
            
        # compute outflow, storage and area
        inflow.name = 'inflow'
        storage = pd.Series(index=inflow.index, dtype=float, name='storage')
        outflow = pd.Series(index=inflow.index, dtype=float, name='outflow')
        for date, I in tqdm(inflow.items()):
            storage[date] = Vo
            Q, V = self.step(
                I, 
                Vo, 
                date.dayofyear,
                P=precipitation[date] if precipitation is not None else None, 
                E=evaporation[date] if evaporation is not None else None,
                D=demand[date] if demand is not None else None
            )
            outflow[date] = Q
            # update current storage
            Vo = V     
                
        return pd.concat((storage, inflow, outflow), axis=1)