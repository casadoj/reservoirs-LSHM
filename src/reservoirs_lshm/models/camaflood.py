import numpy as np
from scipy.special import expit
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from typing import Union, List, Tuple, Dict, Optional
from pathlib import Path
import logging
logger = logging.getLogger(__name__)

from ..utils.plots import reservoir_analysis
from .basemodel import Reservoir

        
class Camaflood(Reservoir):
    """Representation of a reservoir according to Hanazaki, Yamazaki & Yoshimura (2021)."""
    
    def __init__(
        self,
        Vmin: float,
        Vf:float,
        Ve: float,
        Vtot: float,
        Qn: float,
        Qf: float,
        catchment: int,
        Atot: Optional[int] = None,
        timestep: int = 86400
    ):
        """        
        Parameters:
        -----------
        Vmin: float
            Volume (m3) associated to the conservative storage
        Vf: float
            Volume (m3) associated to the flood storage
        Ve: float
            Volume (m3) associated with an emergency situation
        Vtot: float
            Total reservoir storage capacity (m3)
        Qn: float
            Normal outflow (m3/s)
        Qf: float
            Outflow (m3/s) in case of flood
        catchment: integer
            Area (m2) of the reservoir catchment
        Atot: integer (optional)
            Reservoir area (m2) at maximum capacity
        timestep: int
            Simulation time step in seconds.
        """
        
        super().__init__(Vmin, Vtot, Qmin=None, Qf=Qf, Atot=Atot, timestep=timestep)
        
        # storage limits
        self.Vf = Vf
        self.Ve = Ve
        
        # outflow limits
        self.Qn = Qn

        # catchment area
        self.catchment = catchment
        
        # release coefficient
        # self.k = max(1 - (Vtot - Vf) / (catchment * .2), 0)

    def weighing(
        self,
        I: Union[float, np.ndarray],
        Qmin: float,
        eps: float = 5e-3
    ) -> Union[float, np.ndarray]:
        """
        Computes a sigmoid-based weight to smoothly transition between 
        'normal' and 'flood' reservoir release regimes based on inflow.
    
        This function returns a value between ~0 and ~1 that increases 
        smoothly as inflow (I) approaches and exceeds the flood threshold (Qf). 
        The transition replaces sharp mode switching with a smooth blending.
    
        Parameters
        ----------
        I : float or array-like
            Inflow to the reservoir at the current timestep.
        eps : float, optional (default=0.005)
            Accepted error (epsilon). Controls the width of the sigmoid transition. 
            A smaller epsilon makes the transition sharper; a larger factor 
            makes it smoother.
    
        Returns
        -------
        weight : float or ndarray
            A smooth weight in the [0, 1] range that can be used to blend 
            between normal and flood release equations.
        """
        
        # half width of the sigmoid function based on the error
        half_width = np.log((1 - eps) / eps)
        
        # linear transformation: mean and scale
        mu = np.mean([Qmin, self.Qf]).item()
        scale = half_width / (self.Qf - mu)
        x = (I - mu) * scale

        # compute weight: sigmoid
        weight = expit(x)

        return weight
        
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

        # release coefficient
        self.k = max(1 - (self.Vtot - V) / (self.catchment * .2), 0)
        
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
        
        # ouflow depending on the inflow and storage level
        if V < self.Vmin:
            Q = V / self.Vf * self.Qn
        elif V <= self.Ve:
            # Qlow: when inflow is lower than Qf
            # Qhigh: when inflow is greater or equal than Qf
            if V < self.Vf:
                Qlow = self.Vmin / self.Vf * self.Qn + ((V - self.Vmin) / (self.Ve - self.Vmin))**2 * (self.Qf - self.Vmin / self.Vf * self.Qn)
                Qhigh = self.Vmin / self.Vf * self.Qn + (V - self.Vmin) / (self.Vf - self.Vmin) * (self.Qf - self.Vmin / self.Vf * self.Qn)
            elif V < self.Ve:
                Qlow = self.Vmin / self.Vf * self.Qn + ((V - self.Vmin) / (self.Ve - self.Vmin))**2 * (self.Qf - self.Vmin / self.Vf * self.Qn)
                Qhigh = self.Qf + self.k * (V - self.Vf) / (self.Ve - self.Vf) * (I - self.Qf)
            # compute release as a weighed sum of the "normal" and "flood" mode
            w = self.weighing(I, Qlow)
            Q = (1 - w) * Qlow + w * Qhigh
        elif V >= self.Ve:
            Q = self.Qf if I < self.Qf else I
            
            

        # limit outflow so the final storage is between 0 and 1
        eps = 1e-3
        if V - Q * self.timestep > self.Vtot:
            Q = (V - self.Vtot) / self.timestep + eps
        # elif V - Q * self.timestep < self.Vmin:
        #     Q = (V - self.Vmin) / self.timestep - eps if V >= self.Vmin else 0
        elif V - Q * self.timestep < 0:
            Q = V / self.timestep - eps
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
        I: Union[float, pd.Series], 
        modified: bool = True
    ):
        """Given a time series of reservoir storage (m3) and a value or a time series of inflow (m3/s), it computes the ouflow (m3/s). This function is only meant for explanatory purposes; since the volume time series is given, the computed outflow does not update the reservoir storage. If the intention is to simulate the behaviour of the reservoir, refer to the function "simulate"
        
        Parameters:
        -----------
        V: pd.Series
            Time series of reservoir storage (m3)
        I: Union[float, pd.Series]
            Reservor inflow (m3/s)
        modified: bool
            Whether to use the modified (default) of the orinigal Hanazaki's routine. The modified routine avoids the breaks in the outflow function at Vmin, Vf and Ve.
            
        Returns:
        --------
        O: pd.Series
            Time series of reservoir outflow (m3/s)
        """

        if isinstance(I, float) or isinstance(I, int):
            assert I >= 0, '"I" must be a positive value'
            I = pd.Series(I, index=V.index)

        # # release coefficient
        # self.k = max(1 - (self.Vtot - V) / (self.catchment * .2), 0)
        
        maskI = I < self.Qf
        maskV1 = V < self.Vmin
        maskV2 = (self.Vmin <= V) & (V < self.Vf)
        maskV3 = (self.Vf <= V) & (V < self.Ve)
        maskV4 = self.Ve <= V

        O = pd.Series(index=V.index, dtype=float)

        # if in the lower storage level
        O[maskV1] = self.Qn * V[maskV1] / self.Vf

        # if input below flood... 
        # ... and storage below emergency level
        mask = (maskI & maskV2) | (maskI & maskV3)
        if np.sum(mask) > 0:
            if modified:
                O[mask] = self.Vmin / self.Vf * self.Qn + ((V[mask] - self.Vmin) / (self.Ve - self.Vmin))**2 * (self.Qf - self.Vmin / self.Vf * self.Qn)
            else:
                O[mask] = .5 * self.Qn + ((V[mask] - self.Vmin) / (self.Ve - self.Vmin))**2 * (self.Qf - self.Qn)
        # ... and storage above emergency level
        O[maskI & maskV4] = self.Qf

        # if inflow over flood...
        # ... and storage in the normal level
        mask = ~maskI & maskV2
        if np.sum(mask) > 0:
            if modified:
                O[mask] = self.Vmin / self.Vf * self.Qn + (V[mask] - self.Vmin) / (self.Vf - self.Vmin) * (self.Qf - self.Vmin / self.Vf * self.Qn)
            else:
                O[mask] = .5 * self.Qn + (V[mask] - self.Vmin) / (self.Vf - self.Vmin) * (self.Qf - self.Qn)
        # ... and storage in flood level
        mask = ~maskI & maskV3
        if np.sum(mask) > 0:
            O[mask] = self.Qf + self.k * (V[mask] - self.Vf) / (self.Ve - self.Vf) * (I[mask] - self.Qf)
        # ... and storage in emergency level
        O[~maskI & maskV4] = I[~maskI & maskV4]

        return O
    
    def plot_routine(
        self, 
        modified: bool = True, 
        ax: Axes = None, 
        **kwargs
    ):
        """It creates a plot that explains the reservoir routine.
        
        Parameters:
        -----------
        modified: bool
            Whether to use the modified (default) of the orinigal Hanazaki's routine. The modified routine avoids the breaks in the outflow function at Vmin, Vf and Ve.
        ax: Axes
            If provided, the plot will be added to the given axes
        """

        # dummy storage time series
        V = pd.Series(np.linspace(0, self.Vtot + .01, 1000))

        # create scatter plot
        if ax is None:
            fig, ax = plt.subplots(figsize=kwargs.get('figsize', (4, 4)))

        # outflow when inflow is lower than the flood outflow
        outflow1 = self.routine(V, .9 * self.Qf, modified=modified)
        ax.scatter(V, outflow1, s=.05, c='C0', label=r'$I < Q_f$')

        # outflow when inflow is larger than the flood outflow
        outflow2 = self.routine(V, 1.2 * self.Qf, modified=modified)
        ax.scatter(V, outflow2, s=.05, c='C1', label=r'$I \geq Q_f$')

        # reference storages and outflows
        vs = [self.Vmin, self.Vf, self.Ve]
        if modified:
            qs = [self.Vmin / self.Vf * self.Qn, self.Qf, self.Qf]
        else:
            qs = [.5 * self.Qn, self.Qf, self.Qf]
        
        for v, q in zip(vs, qs):
            ax.vlines(v, 0, q, color='k', ls=':', lw=.5, zorder=0)
            ax.hlines(q, 0, v, color='k', ls=':', lw=.5, zorder=0)
        
        # labels
        if modified:
            ax.text(0, qs[0], r'$\frac{V_c}{V_f} Q_n$', ha='left', va='bottom')
        else:
            ax.text(0, qs[0], r'$0.5 Q_n$', ha='left', va='bottom')
        ax.text(0, qs[1], r'$Q_f$', ha='left', va='bottom')
        ax.text(self.Vmin, 0, r'$V_c$', rotation=90, ha='left', va='bottom')
        ax.text(self.Vf, 0, r'$V_f$', rotation=90, ha='right', va='bottom')
        ax.text(self.Ve, 0, r'$V_e$', rotation=90, ha='right', va='bottom')
        
        # setup
        ax.set(xlim=(0, self.Vtot),
               xlabel='storage (hm3)',
               ylim=(0, None),
               ylabel='outflow (m3/s)')
        ax.legend(frameon=False, loc=2)
    
    def get_params(self):
        """It generates a dictionary with the reservoir paramenters in the Hanazaki model."""

        params = {
            'Vmin': self.Vmin,
            'Vf': self.Vf,
            'Ve': self.Ve,
            'Vtot': self.Vtot,
            'Qn': self.Qn,
            'Qf': self.Qf,
            # 'k': self.k,
            'Atot': self.Atot
        }
        # params = {key: float(value) for key, value in params.items()}
        params = {key: float(value) if value is not None else None for key, value in params.items()}

        return params
    
    def compare(
        self, 
        series1: pd.DataFrame, 
        series2: pd.DataFrame = None, 
        save: Union[Path, str] = None, 
        **kwargs
    ):
        """It compares two reservoir timeseries (inflow, outflow and storage) using the function 'reservoir_analysis'. If only 1 time series is given, the plot will simply show the reservoir behaviour of that set of time series.
        
        Inputs:
        -------
        series1: pd.DataFrame
            A table with the time series of 'inflow', 'outflow' and 'storage'
        series2: pd.DataFrame
            A second table with the time series of 'inflow', 'outflow' and 'storage'
        save: Union[Path, str]
            Directory and file where the figure will be saved
            
        kwargs:
        -------
        title: str
            If provided, title of the figure
        labels: List[str]
            A list of 2 strings to be used as labels for each set of time series
        alpha: float
            The transparency of the scatter plot
        """
        
        # storage limits
        Vlims = np.array([self.Vmin, self.Vf, self.Ve]) / self.Vtot
        
        # outflow limits
        Qlims = np.array([self.Vmin / self.Vf * self.Qn, self.Qf, self.Qf]) / self.Qf
        
        # plot analysis
        series1_norm = self.normalize_timeseries(series1)
        if series2 is not None:
            series2_norm = self.normalize_timeseries(series2)
        else:
            series2_norm = series2
        reservoir_analysis(
            series1_norm,
            series2_norm,
            x_thr=Vlims, 
            y_thr=Qlims,
            title=kwargs.get('title', None),
            labels=kwargs.get('labels', ['sim', 'obs']),
            alpha=kwargs.get('alpha', .05),
            save=save
        )
