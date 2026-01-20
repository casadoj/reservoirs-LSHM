import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from tqdm.auto import tqdm
from typing import Union, List, Tuple, Dict, Optional
import logging
logger = logging.getLogger(__name__)

from .basemodel import Reservoir

        
class Lisflood(Reservoir):
    """Representation of a reservoir in the LISFLOOD-OS hydrological model."""
    
    def __init__(
        self,
        Vmin: float,
        Vn: float,
        Vn_adj: float,
        Vf: float,
        Vtot: float,
        Qmin: float,
        Qn: float,
        Qf: float,
        k: float = 1.2,
        Atot: Optional[int] = None,
        timestep: int = 86400
    ):
        """
        Parameters:
        -----------
        Vmin: float
            Volume (m3) associated to the conservative storage
        Vn: float
            Volume (m3) associated to the normal storage
        Vn_adj: float
            Volume (m3) associated to the adjusted (calibrated) normal storage
        Vf: float
            Volume (m3) associated to the flood storage
        Vtot: float
            Total reservoir storage capacity (m3)
        Qmin: float
            Minimum outflow (m3/s)
        Qn: float
            Normal outflow (m3/s)
        Qf: float
            Non-damaging outflow (m3/s)
        k: float
            Flood release factor. It allows for releases k times the inflow, i.e., larger than the inflow
        Atot: integer (optional)
            Reservoir area (m2) at maximum capacity
        timestep: int
            Simulation time step in seconds.
        """
        
        super().__init__(Vmin, Vtot, Qmin, Qf, Atot, timestep)
        
        # storage limits
        self.Vn = Vn
        self.Vn_adj = Vn_adj
        self.Vf = Vf
        
        # outflow limits
        self.Qn = Qn
        self.k = k
    
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
            V = np.max([0, V - E * 1e-3 * A])
        if D:
            # demand can't withdraw water below the minimum storage
            V = np.max([self.Vmin, V - D])
        
        # ouflow depending on the storage level
        if V < 2 * self.Vmin:
            Q = self.Qmin
        elif V < self.Vn:
            Q = self.Qmin + (self.Qn - self.Qmin) * (V - 2 * self.Vmin) / (self.Vn - 2 * self.Vmin)
        elif V < self.Vn_adj:
            Q = self.Qn
        elif V < self.Vf:
            Q = self.Qn + (self.Qf - self.Qn) * (V - self.Vn_adj) / (self.Vf - self.Vn_adj)
            if Q > self.k * I:
                Q = np.max([self.k * I, self.Qn])
        elif V > self.Vf:
            Q = np.max([(V - self.Vf) / self.timestep, np.min([self.Qf, np.max([self.k * I, self.Qn])])])
        
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
        
    def routine(
        self,
        V: pd.Series,
        I: Union[float, pd.Series]
    ) -> pd.Series:
        """Given a time series of reservoir storage (m3) and a value or a time series of inflow (m3/s), it computes the ouflow (m3/s). This function is only meant for explanatory purposes; since the volume time series is given, the computed outflow does not update the reservoir storage. If the intention is to simulate the behaviour of the reservoir, refer to the function "simulate"
        
        Parameters:
        -----------
        V: pd.Series
            Time series of reservoir storage (m3)
        I: Union[float, pd.Series]
            Reservor inflow (m3/s)
            
        Returns:
        --------
        O: pd.Series
            Time series of reservoir outflow (m3/s)
        """
        
        if isinstance(I, float) or isinstance(I, int):
            assert I >= 0, '"I" must be a positive value'
            I = pd.Series(I, index=V.index)
        
        O1 = V / self.timestep 
        O1[O1 > self.Qmin] = self.Qmin
        O = O1.copy()
        
        O2 = self.Qmin + (self.Qn - self.Qmin) * (V - 2 * self.Vmin) / (self.Vn - 2 * self.Vmin)
        maskV2 = (2 * self.Vmin <= V) & (V < self.Vn)
        O[maskV2] = O2[maskV2]
        
        O3 = pd.Series(self.Qn, index=V.index)
        maskV3 = (self.Vn <= V) & (V < self.Vn_adj)
        O[maskV3] = O3[maskV3]
        
        O4 = self.Qn + (self.Qf - self.Qn) * (V - self.Vn_adj) / (self.Vf - self.Vn_adj)
        maskV4 = (self.Vn_adj <= V) & (V < self.Vf)
        O[maskV4] = O4[maskV4]
        
        Omax = 1.2 * I
        Omax[Omax < self.Qn] = self.Qn
        Omax[Omax > self.Qf] = self.Qf
        O5 = pd.concat(((V - self.Vf - .01 * self.Vtot) / self.timestep, Omax), axis=1).max(axis=1)
        maskV5 = self.Vf <= V
        O[maskV5] = O5[maskV5]
        
        Oreg = I
        Oreg[Oreg < self.Qn] = self.Qn
        Oreg = pd.concat((O, Oreg), axis=1).min(axis=1)
        maskO = (O > 1.2 * I) & (O > self.Qn) & (V < self.Vf)
        O[maskO] = Oreg[maskO]
        
        temp = pd.concat((O1, O2, O3, O4, O5, Omax, Oreg), axis=1)
        temp.columns = ['O1', 'O2', 'O3', 'O4', 'O5', 'Omax', 'Oreg']
        self.O = temp
        
        return O
       
    def plot_routine(
        self, 
        ax: Axes = None, 
        **kwargs
    ):
        """It creates a plot that explains the reservoir routine.
        
        Parameters:
        -----------
        ax: Axes
            If provided, the plot will be added to the given axes
        """

        # dummy storage time series
        V = pd.Series(np.linspace(0, self.Vtot + .01, 1000))

        # create scatter plot
        if ax is None:
            fig, ax = plt.subplots(figsize=kwargs.get('figsize', (5, 5)))

        # outflow
        outflow = self.routine(V, I=self.Qf)
        ax.plot(V, outflow, lw=1, c='C0')

        # reference storages and outflows
        vs = [self.Vmin, 2 * self.Vmin, self.Vn, self.Vn_adj, self.Vf]
        qs = [self.Qmin, self.Qmin, self.Qn, self.Qn, self.Qf]
        for v, q in zip(vs, qs):
            ax.vlines(v, 0, q, color='k', ls=':', lw=.5, zorder=0)
            ax.hlines(q, 0, v, color='k', ls=':', lw=.5, zorder=0)
        
        # labels
        ax.text(0, self.Qmin, r'$Q_{min}$', ha='left', va='bottom')
        ax.text(0, self.Qn, r'$Q_{n,adj}$', ha='left', va='bottom')
        ax.text(0, self.Qf, r'$Q_nd$', ha='left', va='bottom')
        ax.text(self.Vn, 0, r'$V_n$', rotation=90, ha='right', va='bottom')
        ax.text(self.Vn_adj, 0, r'$V_{n,adj}$', rotation=90, ha='right', va='bottom')
        ax.text(self.Vf, 0, r'$V_f$', rotation=90, ha='right', va='bottom')
        
        # setup
        ax.set(xlim=(0, self.Vtot),
               xlabel='storage (hm3)',
               ylim=(0, None),
               ylabel='outflow (m3/s)')
        ax.set_title('LISFLOOD reservoir routine')
        
    def get_params(self) -> Dict:
        """It generates a dictionary with the reservoir parameters
        
        Returns:
        --------
        params: Dict
            A dictionary with the name and value of the reservoir parameters
        """

        params = {
            'Vmin': self.Vmin,
            'Vn': self.Vn,
            'Vn_adj': self.Vn_adj,
            'Vf': self.Vf,
            'Vtot': self.Vtot,
            'Qmin': self.Qmin,
            'Qn': self.Qn,
            'Qf': self.Qf,
            'k': self.k,
            'Atot': self.Atot
        }
        # params = {key: float(value) for key, value in params.items()}
        params = {key: float(value) if value is not None else None for key, value in params.items()}

        return params
    

def plot_lisflood(reservoir, timeseries, ax=None, **kwargs):

    lw = kwargs.get('lw', 1.2)
    c = kwargs.get('c', 'steelblue')
    figsize = kwargs.get('figsize', (4, 4))

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    
    V = np.array([0, 2 * reservoir.Vmin])
    Q = np.zeros_like(V) + reservoir.Qmin
    ax.plot(V * 1e-6, Q, lw=lw, c='steelblue')
    
    V = np.linspace(2 * reservoir.Vmin, reservoir.Vn, 100)
    Q = reservoir.Qmin + (reservoir.Qn - reservoir.Qmin) * (V - 2 * reservoir.Vmin) / (reservoir.Vn - 2 * reservoir.Vmin)
    ax.plot(V * 1e-6, Q, lw=lw, c=c)
    
    V = np.linspace(reservoir.Vn, reservoir.Vn_adj, 100)
    Q = np.zeros_like(V) + reservoir.Qn
    ax.plot(V * 1e-6, Q, lw=lw, c=c)
    
    V = np.linspace(reservoir.Vn_adj, reservoir.Vf, 100)
    Q = reservoir.Qn + (reservoir.Qf - reservoir.Qn) * (V - reservoir.Vn_adj) / (reservoir.Vf - reservoir.Vn_adj)
    #if Q > reservoir.k * I:
    #    Q = np.max([reservoir.k * I, reservoir.Qn])
    ax.plot(V * 1e-6, Q, lw=lw, c=c)
    
    #V = np.linspace(reservoir.Vf, reservoir.Vtot, 100)
    #I = 1.1 * reservoir.Qf
    #Q = np.maximum((V - reservoir.Vf) / reservoir.timestep, np.min([reservoir.Qf, np.max([reservoir.k * I, reservoir.Qn])]))
    #ax.plot(V * 1e-6, Q, lw=lw, c=c)
    
    ax.scatter(timeseries.storage * 1e-6, timeseries.outflow, marker='.', s=.5, color='lightgrey', alpha=.5, zorder=0)
    
    # reference storages and outflows
    vs = [2 * reservoir.Vmin, reservoir.Vn, reservoir.Vn_adj, reservoir.Vf]
    qs = [reservoir.Qmin, reservoir.Qn, reservoir.Qn, reservoir.Qf]
    for v, q in zip(vs, qs):
        ax.vlines(v * 1e-6, 0, q, color='k', ls=':', lw=.5, zorder=0)
        ax.hlines(q, 0, v * 1e-6, color='k', ls=':', lw=.5, zorder=0)
                
    # labels
    ax.text(0, reservoir.Qmin, r'$Q_c$', ha='left', va='bottom')
    ax.text(0, reservoir.Qn, r'$Q_n$', ha='left', va='bottom')
    ax.text(0, reservoir.Qf, r'$Q_f$', ha='left', va='bottom')
    ax.text(reservoir.Vn * 1e-6, 0, r'$V_n$', rotation=90, ha='right', va='bottom')
    ax.text(reservoir.Vn_adj * 1e-6, 0, r"$V_n'$", rotation=90, ha='right', va='bottom')
    ax.text(reservoir.Vf * 1e-6, 0, r'$V_f$', rotation=90, ha='right', va='bottom')
    
    # setup
    ax.set(xlim=(0, reservoir.Vtot * 1e-6),
           xlabel=r'Storage [$hm^3$]',
           ylim=(0, None),
           ylabel=r'Outflow [$m^3/s$]')
    ax.set_title('LISFLOOD')
