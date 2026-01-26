#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import logging
import yaml
import pickle
from pathlib import Path
import pandas as pd
from tqdm import tqdm
import spotpy
from datetime import datetime
import sys # Added for sys.platform
import os # Added for setting environment variables if needed

# Import MPI
from mpi4py import MPI

# Assuming these are available from your project structure
# from . import Config, read_attributes, read_timeseries
# from .models import get_model
# from .calibration import get_calibrator, read_results
# from .utils.metrics import KGEmod, compute_performance
# from .utils.timeseries import create_demand
# from .utils.logging import setup_logger

# Mock imports for standalone example
class Config:
    def __init__(self, config_path):
        self.MODEL = 'camaflood'
        self.PATH_CALIB = Path('./calibration_results')
        self.RESERVOIRS_FILE = 'reservoirs.csv'
        self.PATH_DATA = Path('./data')
        self.PERIODS_FILE = 'periods.pkl'
        self.INFLOW = 'inflow'
        self.PRECIPITATION = 'precipitation'
        self.EVAPORATION = 'evaporation'
        self.DEMAND = 'demand'
        self.COMPLEXES = 4 # This is the ngs value for SPOTPY, not the total MPI processes
        self.MAX_ITER = 100
        self.SPINUP = 0
        self.TARGET = 'outflow'
        self.PATH_DEF = Path('./default_sims')
        self.PATH_CALIB.mkdir(parents=True, exist_ok=True)
        self.PATH_DEF.mkdir(parents=True, exist_ok=True) # Ensure default_sims directory exists for mock

def setup_logger(name, log_level, log_file):
    logging.basicConfig(level=log_level, filename=log_file, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    return logging.getLogger(name)

# Mock data and functions for demonstration
def read_attributes(path, reservoirs):
    return pd.DataFrame({'CATCH_SKM': {r: 1000 for r in reservoirs}, 'AREA_SKM': {r: 100 for r in reservoirs}}, index=reservoirs)

def read_timeseries(path, reservoirs, periods, variables):
    # Create dummy timeseries data for demonstration
    timeseries_data = {}
    date_range = pd.date_range(start='2000-01-01', periods=100, freq='D')
    for r_id in reservoirs:
        df = pd.DataFrame(index=date_range)
        for var in variables:
            df[var] = 10 + (r_id * 0.1) + pd.np.random.rand(len(date_range)) # Use pd.np for numpy in pandas context
        df['outflow'] = df['inflow'] * 0.8 # Simple relationship
        df['storage'] = df['outflow'].cumsum() # Simple relationship for storage
        timeseries_data[r_id] = df
    return timeseries_data

def get_calibrator(model, **kwargs):
    # Mock calibrator
    class MockCalibrator:
        def __init__(self):
            self.parameters = {'param1': 0.5, 'param2': 0.1} # Example parameters

        def pars2attrs(self, params):
            return {'param1': params[0], 'param2': params[1]}

        def simulation(self, evaluation, parameters):
            # Simulate some data based on parameters
            simulated_outflow = evaluation * parameters['param1'] # Simplified
            return simulated_outflow.values # SPOTPY expects numpy array

        def parameters(self):
            # SPOTPY requires this to define parameters
            return spotpy.parameter.Uniform('param1', 0, 1), spotpy.parameter.Uniform('param2', 0, 0.5)

    return MockCalibrator()

def read_results(filename):
    # Mock read_results
    return None, {'param1': 0.6, 'param2': 0.2}

def get_model(model, **kwargs):
    # Mock model
    class MockModel:
        def get_params(self):
            return {'param1': 0.6, 'param2': 0.2}

        def simulate(self, inflow, Vo, precipitation, evaporation, demand):
            # Mock simulation, return a DataFrame
            date_range = pd.date_range(start='2000-01-01', periods=100, freq='D')
            return pd.DataFrame({'outflow': inflow * 0.7}, index=date_range)

        def scatter(self, sim, ts, **kwargs):
            # Mock scatter plot
            print(f"Generating scatter plot for {kwargs.get('title')}")

        def lineplot(self, sim, ts, **kwargs):
            # Mock line plot
            print(f"Generating line plot for {kwargs.get('save')}")
    return MockModel()

def KGEmod(simulated, observed):
    # Mock KGEmod
    return 0.85

def compute_performance(obs, sim):
    # Mock performance computation
    return pd.Series({'KGEmod': 0.85, 'RMSE': 10.0})

def create_demand(outflow, water_stress, window):
    # Mock create_demand
    return outflow * water_stress

# End of mock imports

def main():
    # Initialize MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # parse arguments (only on rank 0, then broadcast if needed, or parse independently)
    parser = argparse.ArgumentParser(
        description="""
        Run the calibration script with a specified configuration file.
        It calibrates the reservoir model parameters of the defined routine using the
        SCE-UA (Shuffle Complex Evolution-University of Arizona) algorithm for each of
        the selected reservoirs.
        The optimal parameters are simulated and plotted, if possible comparing against
        a simulation with default parameters
        """
    )
    parser.add_argument('-c', '--config-file', type=str, required=True, help='Path to the YAML configuration file (e.g., config.yml).')
    parser.add_argument('-o', '--overwrite', action='store_true', default=False, help='Overwrite existing calibration results.')
    args = parser.parse_args()

    # read configuration file
    cfg = Config(args.config_file) # Each process reads its own config

    # set up logger (each process logs to a unique file to avoid conflicts)
    log_file_name = f'{datetime.now():%Y%m%d%H%M%S}_calibrate_{cfg.MODEL}_rank{rank}.log'
    logger = setup_logger(
        name=__name__,
        log_level=logging.INFO,
        log_file=log_file_name
    )

    logger.info(f'Calibration results will be saved in: {cfg.PATH_CALIB}')
    logger.info(f'MPI Rank {rank}/{size} started.')

    # === Load reservoir list ===
    # Only rank 0 loads the full list, then distributes
    if rank == 0:
        try:
            reservoirs_full_list = pd.read_csv(cfg.RESERVOIRS_FILE, header=None).squeeze().tolist()
            logger.info(f'Full list of reservoirs: {len(reservoirs_full_list)}')
        except IOError:
            logger.exception(f'Failed to open {cfg.RESERVOIRS_FILE}')
            comm.Abort() # Abort all processes on error
            raise
    else:
        reservoirs_full_list = None

    # Broadcast the full list of reservoirs to all processes
    reservoirs_full_list = comm.bcast(reservoirs_full_list, root=0)

    # Distribute reservoirs among processes
    # Simple chunking for demonstration
    # You might want more sophisticated dynamic load balancing for real applications
    my_reservoirs = []
    for i, res_id in enumerate(reservoirs_full_list):
        if i % size == rank:
            my_reservoirs.append(res_id)
    
    logger.info(f'Rank {rank} assigned {len(my_reservoirs)} reservoirs: {my_reservoirs}')

    # === Load attributes ===
    # Each process reads attributes for all reservoirs that *might* be relevant,
    # or just for its assigned reservoirs if attributes are large.
    # For simplicity, here we read all if needed.
    try:
        attributes = read_attributes(cfg.PATH_DATA / 'attributes', reservoirs_full_list)
        logger.info(f'{attributes.shape[0]} reservoirs in the attribute tables on Rank {rank}')
    except IOError:
        logger.exception('Failed to read attribute tables from {0} on Rank {1}'.format(cfg.PATH_DATA / 'attributes', rank))
        comm.Abort()
        raise

    # === Load time periods ===
    try:
        with open(cfg.PERIODS_FILE, 'rb') as file:
            periods = pickle.load(file)
    except IOError:
        logger.exception(f'Failed to open {cfg.PERIODS_FILE} on Rank {rank}')
        comm.Abort()
        raise

    # === read time series ===
    try:
        inputs = [var for var in [cfg.INFLOW, cfg.PRECIPITATION, cfg.EVAPORATION, cfg.DEMAND] if var]
        outputs = ['storage', 'outflow']
        timeseries = read_timeseries(
            path=cfg.PATH_DATA / 'time_series' / 'csv',
            reservoirs=reservoirs_full_list, # Read all timeseries
            periods=periods,
            variables=inputs + outputs
        )
        logger.info(f'{len(timeseries)} reservoirs with timeseries on Rank {rank}')
    except IOError:
        logger.exception('Failed to read time series from {0} on Rank {1}'.format(cfg.PATH_DATA / 'time_series' / 'csv', rank))
        comm.Abort()
        raise

    # === Calibration ===
    # Iterate only over reservoirs assigned to this rank
    for grand_id in tqdm(my_reservoirs, desc=f'Rank {rank} calibrating reservoir'):

        # Ensure output directories are unique for each process or handle file naming carefully
        # In this case, dbname already includes grand_id, so it should be fine.
        dbname = f'{cfg.PATH_CALIB}/{grand_id}_samples'
        if Path(f'{dbname}.csv').is_file() and not args.overwrite:
            logger.info(f'Rank {rank}: Calibration already exists for reservoir {grand_id}, skipping (use --overwrite to force)')
            continue
        
        logger.info(f'Rank {rank}: Calibrating reservoir {grand_id}')
        
        ts = timeseries[grand_id] # Get timeseries for the assigned grand_id
        
        # define input time series
        inflow = ts[cfg.INFLOW]
        precipitation = ts[cfg.PRECIPITATION] if cfg.PRECIPITATION in ts.columns else None
        evaporation = ts[cfg.EVAPORATION] if cfg.EVAPORATION in ts.columns else None
        demand = ts[cfg.DEMAND] if cfg.DEMAND in ts.columns else None
        if cfg.MODEL == 'mhm':
            bias = ts.outflow.mean() / inflow.mean()
            demand = create_demand(
                ts.outflow,
                water_stress=min(1, bias),
                window=28
            )
            
        # storage attributes (m3)
        Vtot = ts.storage.max()
        Vmin = max(0, min(0.1 * Vtot, ts.storage.min()))
        # flow attributes (m3/s)
        Qmin = max(0, ts.outflow.min())
        # catchment area (m2)
        catchment = int(attributes.loc[grand_id, 'CATCH_SKM'] * 1e6) if cfg.MODEL == 'camaflood' else None
        # reservoir area (m2)
        Atot = int(attributes.loc[grand_id, 'AREA_SKM'] * 1e6)

        # calibrate
        try:
            # configure calibration kwargs
            cal_cfg = {}
            if cfg.MODEL == 'camaflood':
                cal_cfg.update({'catchment': catchment})
            # elif cfg.MODEL == 'mhm':
            #    cal_cfg.update({'demand': demand})
            
            # initialize the calibration setup
            calibrator = get_calibrator(
                cfg.MODEL,
                inflow=inflow,
                storage=ts.storage,
                outflow=ts.outflow,
                precipitation=precipitation,
                evaporation=evaporation,
                demand=demand,
                Vmin=Vmin,
                Vtot=Vtot,
                Qmin=Qmin,
                Atot=Atot,
                target=cfg.TARGET,
                obj_func=KGEmod,
                spinup=cfg.SPINUP,
                **cal_cfg
            )
            
            # SPOTPY expects a unique `parallel='mpi'` instance per rank.
            # If your `cfg.COMPLEXES` is meant for SPOTPY's *internal* parallelism (e.g., within one SCE-UA run),
            # then keep `parallel='mpi'` here. If `cfg.COMPLEXES` is supposed to be the *total* MPI processes
            # for your script, and each rank should run SCE-UA sequentially, then remove `parallel='mpi'`.
            # Based on the problem description, you have an outer MPI loop and SPOTPY's inner MPI, which is complex.
            # For simplicity, let's assume `cfg.COMPLEXES` is for SPOTPY's internal parallelization,
            # and our `mpi4py` handles the outer parallelization.
            # If `cfg.COMPLEXES` should be `size` for `spotpy`, then set `ngs=size` below.
            
            # IMPORTANT: If cfg.COMPLEXES is intended for the _total_ MPI processes for the entire script,
            # and each process is doing an independent calibration, then you should NOT use parallel='mpi' here
            # for SPOTPY, or you should set ngs = 1 for SPOTPY.
            # Assuming cfg.COMPLEXES is still for SPOTPY's internal complex number
            sceua = spotpy.algorithms.sceua(
                calibrator,
                dbname=dbname,
                dbformat='csv',
                parallel='mpi', # This tells SPOTPY to use MPI for *its* internal calculations
                save_sim=False,
                # seed=42
            )
            # launch calibration
            sceua.sample(
                cfg.MAX_ITER,
                ngs=cfg.COMPLEXES, # This is the number of complexes for SCE-UA
                kstop=3,
                pcento=0.01,
                peps=0.1
            )
            logger.info(f'Rank {rank}: Calibration of reservoir {grand_id} successfully finished')
        except RuntimeError:
            logger.exception(f'Rank {rank}: Reservoir {grand_id} could not be calibrated')
            continue
            
        # simulate optimized reservoir
        try:
            # read calibration results
            results, parameters = read_results(f'{dbname}.csv')
            
            # convert parameter into reservoir attributes
            calibrated_attrs = calibrator.pars2attrs(list(parameters.values()))

            # declare the reservoir with optimal parameters
            res = get_model(cfg.MODEL, **calibrated_attrs)

            # export calibrated parameters
            with open(cfg.PATH_CALIB / f'{grand_id}_optimal_parameters.yml', 'w') as file:
                yaml.dump(res.get_params(), file)

            # simulate the reservoir
            Vo = ts.storage.iloc[0]
            sim_cal = res.simulate(
                inflow=inflow,
                Vo=None if pd.isna(Vo) else Vo,
                precipitation=precipitation,
                evaporation=evaporation,
                demand=demand,
            )
            sim_cal.to_csv(cfg.PATH_CALIB / f'{grand_id}_simulation.csv', float_format='%.3f')
            
            logger.info(f'Rank {rank}: Simulation of the calibrated reservoir {grand_id} successfully finished')
            
        except RuntimeError:
            logger.exception(f'Rank {rank}: Calibrated reservoir {grand_id} could not be simulated')
            continue
            
        # === Analyse results ===
        
        # performance
        try:
            performance_cal = compute_performance(ts.iloc[cfg.SPINUP:], sim_cal.iloc[cfg.SPINUP:])
            performance_cal.to_csv(cfg.PATH_CALIB / f'{grand_id}_performance.csv', float_format='%.3f')
            logger.info(f'Rank {rank}: Performance of reservoir {grand_id} has been computed')
        except IOError:
            logger.exception(f'Rank {rank}: The performance of reservoir {grand_id} could not be exported')
            
        # scatter plot calibration vs observation
        try:
            res.scatter(
                sim_cal,
                ts,
                spinup=cfg.SPINUP,
                norm=False,
                title=f'grand_id: {grand_id}',
                save=cfg.PATH_CALIB / f'{grand_id}_scatter.jpg',
            )
            logger.info(f'Rank {rank}: Scatter plot of simulation from reservoir {grand_id}')
        except IOError:
            logger.exception(f'Rank {rank}: The scatter plot of reservoir {grand_id} could not be generated')
            
        # line plot calibration (vs default simulation) vs observation
        try:
            file_default = cfg.PATH_DEF / f'{grand_id}_simulation.csv'
            if file_default.is_file():
                sim_def = pd.read_csv(cfg.PATH_DEF / f'{grand_id}_simulation.csv', parse_dates=True, index_col=0)
                sim = {
                    'default': sim_def,
                    'calibrated': sim_cal
                }
            else:
                sim = {'calibrated': sim_cal}
            res.lineplot(
                sim,
                ts,
                spinup=cfg.SPINUP,
                figsize=(12, 6),
                save=cfg.PATH_CALIB / f'{grand_id}_line.jpg',
            )
            logger.info(f'Rank {rank}: Line plot of simulation from reservoir {grand_id}')
        except IOError:
            logger.exception(f'Rank {rank}: The line plot of reservoir {grand_id} could not be generated')
            
        del res, calibrator, sim_cal, calibrated_attrs, performance_cal #, sim_cfg
        try:
            del sceua
        except NameError: # Use NameError if sceua might not be defined due to earlier error
            pass

    # Synchronize all processes before exiting
    comm.Barrier()
    if rank == 0:
        logger.info("All calibration tasks completed across all ranks.")

if __name__ == "__main__":
    main()