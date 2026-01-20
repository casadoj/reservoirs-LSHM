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

from . import Config, read_attributes, read_timeseries
from .models import get_model
from .calibration import get_calibrator, read_results
from .utils.metrics import KGEmod, compute_performance
from .utils.timeseries import create_demand
from .utils.logging import setup_logger


def main():

    # parse arguments
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
    # parser.add_argument('-p', '--parallel', action='store_true', default=False, help='Parallelize calibration using MPI')
    args = parser.parse_args()

    # read configuration file
    cfg = Config(args.config_file)

    # set up logger
    logger = setup_logger(
        name=__name__,
        log_level=logging.INFO,
        log_file=f'{datetime.now():%Y%m%d%H%M%S}_calibrate_{cfg.MODEL}.log'
    )
    
    logger.info(f'Calibration results will be saved in: {cfg.PATH_CALIB}')
    
    # === Load reservoir list ===
    try:
        reservoirs = pd.read_csv(cfg.RESERVOIRS_FILE, header=None).squeeze().tolist()
    except IOError:
        logger.exception(f'Failed to open {cfg.RESERVOIRS_FILE}')
        raise

    # === Load attributes ===
    try:
        attributes = read_attributes(cfg.PATH_DATA / 'attributes', reservoirs)
        logger.info(f'{attributes.shape[0]} reservoirs in the attribute tables')
    except IOError:
        logger.exception('Failed to read attribute tables from {0}'.format(cfg.PATH_DATA / 'attributes'))
        raise

    # === Load time periods ===
    try:
        with open(cfg.PERIODS_FILE, 'rb') as file:
            periods = pickle.load(file)
    except IOError:
        logger.exception(f'Failed to open {cfg.PERIODS_FILE}')
        raise

    # === read time series ===
    try:
        inputs = [var for var in [cfg.INFLOW, cfg.PRECIPITATION, cfg.EVAPORATION, cfg.DEMAND] if var]
        outputs = ['storage', 'outflow']
        timeseries = read_timeseries(
            path=cfg.PATH_DATA / 'time_series' / 'csv',
            reservoirs=attributes.index,
            periods=periods,
            variables=inputs + outputs
        )
        logger.info(f'{len(timeseries)} reservoirs with timeseries')
    except IOError:
        logger.exception('Failed to read time series from {0}'.format(cfg.PATH_DATA / 'time_series' / 'csv'))
        raise

    # === Calibration ===
    for grand_id, ts in tqdm(timeseries.items(), desc='simulating reservoir'):
            
        dbname = f'{cfg.PATH_CALIB}/{grand_id}_samples'
        if Path(f'{dbname}.csv').is_file() and not args.overwrite:
            logger.info(f'Calibration already exists for reservoir {grand_id}, skipping (use --overwrite to force)')
            continue
        
        logger.info(f'Calibrating reservoir {grand_id}')
        
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
        Vtot = max(attributes.loc[grand_id, 'CAP_MCM'].item() * 1e6, ts.storage.max())
        # Vtot = ts.storage.max()
        Vmin = max(0, min(0.1 * Vtot, ts.storage.min()))
        # flow attributes (m3/s)
        Qmin = max(0, ts.outflow.min())
        # catchment area (m2)
        catchment = int(attributes.loc[grand_id, 'CATCH_SKM'].item() * 1e6) if cfg.MODEL == 'camaflood' else None
        # reservoir area (m2)
        Atot = int(attributes.loc[grand_id, 'AREA_SKM'].item() * 1e6)

        # calibrate
        try:
            # configure calibration kwargs
            cal_cfg = {}
            if cfg.MODEL == 'camaflood':
                cal_cfg.update({'catchment': catchment})
            # elif cfg.MODEL == 'mhm':
            #     cal_cfg.update({'demand': demand})
            # initialize the calibration setup
            calibrator = get_calibrator(
                cfg.MODEL,
                parameters=cfg.PARAMETERS,
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
            # define the sampling method
            sceua = spotpy.algorithms.sceua(
                calibrator, 
                dbname=dbname, 
                dbformat='csv', 
                # parallel='mpi' if args.parallel else 'seq',
                save_sim=False,
                # seed=42
            )
            # launch calibration
            sceua.sample(
                cfg.MAX_ITER, 
                ngs=cfg.COMPLEXES, 
                kstop=cfg.KSTOP, 
                pcento=cfg.PCENTO, 
                peps=cfg.PEPS,
            )
            logger.info(f'Calibration of reservoir {grand_id} successfully finished')
        except RuntimeError:
            logger.exception(f'Reservoir {grand_id} could not be calibrated')
            continue
            
        # simulate optimized reservoir
        try:
            
            # read calibration results
            results, parameters = read_results(f'{dbname}.csv')
            
            # convert parameter into reservoir attributes
            calibrated_attrs = calibrator.pars2attrs(list(parameters.values()))

            # declare the reservoir with optimal parameters
            res = get_model(cfg.MODEL, **calibrated_attrs)
            # if cfg.MODEL == 'camaflood':
            #     res.k = parameters['k']

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
            
            logger.info(f'Simulation of the calibrated reservoir {grand_id} successfully finished')
            
        except RuntimeError:
            logger.exception(f'Calibrated reservoir {grand_id} could not be simulated')
            continue
            
        # === Analyse results ===
        
        # performance
        try:
            performance_cal = compute_performance(ts.iloc[cfg.SPINUP:], sim_cal.iloc[cfg.SPINUP:])
            performance_cal.to_csv(cfg.PATH_CALIB / f'{grand_id}_performance.csv', float_format='%.3f')
            logger.info(f'Performance of reservoir {grand_id} has been computed')
        except IOError:
            logger.exception(f'The performance of reservoir {grand_id} could not be exported')
            
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
            logger.info(f'Scatter plot of simulation from reservoir {grand_id}')
        except IOError:
            logger.exception(f'The scatter plot of reservoir {grand_id} could not be generated')
            
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
            logger.info(f'Line plot of simulation from reservoir {grand_id}')
        except IOError:
            logger.exception(f'The line plot of reservoir {grand_id} could not be generated')
            
        del res, calibrator, sim_cal, calibrated_attrs, performance_cal#, sim_cfg
        try:
            del sceua
        except:
            pass

if __name__ == "__main__":
    main()
