#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import logging
import pickle
from pathlib import Path

import pandas as pd
from tqdm import tqdm
from datetime import datetime

from . import Config, read_attributes, read_timeseries
from .models.starfit.starfit import Starfit
from .utils.metrics import compute_performance
from .utils.logging import setup_logger


def main():

    # parse arguments
    parser = argparse.ArgumentParser(
        description="""
        Run Starfit simulation with the paremeter fitted using `fit_starfit`.
        The simulated time series are saved as CSV files. To analyse the results,
        the code creates a CSV file of performance metrics, and a scatter and a 
        line plot comparing the observed and simulated time series.
        """
    )
    parser.add_argument('-c', '--config-file', type=str, required=True, help='Path to the configuration file')
    parser.add_argument('-o', '--overwrite', action='store_true', help='Overwrite existing simulation files', default=False)
    args = parser.parse_args()

    # set up logger
    logger = setup_logger(
        name=__name__,
        log_level=logging.INFO,
        log_file=f'{datetime.now():%Y%m%d%H%M%S}_run_starfit.log'
    )

    # read configuration file
    cfg = Config(args.config_file)
    logger.info(f'Results will be saved in: {cfg.PATH_CALIB}')
    
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
        logger.exception(f'Failed to read attribute tables')
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
        logger.exception(f'Failed to read time series')
        raise

    # === Simulation ===
    for grand_id, obs in tqdm(timeseries.items(), desc='reservoirs'):

        out_file = cfg.PATH_CALIB / f'{grand_id}_simulation.csv'
        if out_file.exists() and not args.overwrite:
            logger.info(f'Simulation already exists for {grand_id}, skipping (use --overwrite to force)')
            continue
        
        logger.info(f'Simulating reservoir {grand_id}')
        # add day of the year to the observed time series
        obs['doy'] = obs.index.dayofyear

        # load fitted storage model
        try:
            with open(cfg.PATH_DEF.parent / 'NOR' / f'{grand_id}.pkl', 'rb') as file:
                model_storage = pickle.load(file)
            NOR = pd.DataFrame({
                'flood': model_storage["NOR upper bound"],
                'conservation': model_storage["NOR lower bound"]
            })
            Vtot = model_storage['capacity (MCM)'] * 1e6
            logger.debug('Storage model correctly loaded')
        except Exception:
            logger.exception(f'Storage model could not be loaded for {grand_id}')
            continue

        # load fitted release model
        try:
            with open(cfg.PATH_DEF.parent / 'release' / f'{grand_id}.pkl', 'rb') as file:
                model_release = pickle.load(file)
            avg_inflow = model_release['mean inflow (MCM/wk)'] * 1e6 / (7 * 86400)
            Qmin, Qmax = (model_release['constraints'] + 1) * avg_inflow
            logger.debug('Release model correctly loaded')
        except Exception:
            logger.exception(f'Release model could not be loaded for {grand_id}')
            continue

        # declare reservoir
        try:
            res = Starfit(
                Vtot=Vtot,
                avg_inflow=avg_inflow,
                pars_Vf=NOR['flood'],
                pars_Vc=NOR['conservation'],
                pars_Qharm=model_release['harmonic parameters'],
                pars_Qresid=model_release['residual parameters'],
                Qmin=Qmin,
                Qmax=Qmax
            )
            logger.debug('Starfit class correctly declared')
        except Exception:
            logger.exception(f'Could not initialize Starfit for {grand_id}')
            continue

        # run simulation
        try:
            sim = res.simulate(obs.inflow, obs.storage.iloc[0])
            sim.to_csv(out_file, float_format='%.3f')
            logger.info(f'Simulation completed for {grand_id}')
        except Exception:
            logger.exception(f'Simulation failed for {grand_id}')
            continue

        # performance
        try:
            performance_cal = compute_performance(obs, sim)
            performance_cal.to_csv(cfg.PATH_CALIB / f'{grand_id}_performance.csv', float_format='%.3f')
            logger.debug(f'Performance of reservoir {grand_id} computed')
        except Exception:
            logger.exception(f'Could not save performance for {grand_id}')

        # scatter plot simulation vs observation
        try:
            res.scatter(
                sim,
                obs,
                norm=False,
                title=f'grand_id: {grand_id}',
                save=cfg.PATH_CALIB / f'{grand_id}_scatter.jpg',
            )
            logger.debug('Scatter plot successfully generated')
        except Exception:
            logger.exception(f'Failed to generate scatter plot for {grand_id}')

        try:
            res.lineplot(
                {'starfit': sim},
                obs,
                Vlims=[res.Vtot],
                Qlims=[res.Qmin, res.Qmax],
                figsize=(12, 6),
                save=cfg.PATH_CALIB / f'{grand_id}_line.jpg',
            )
            logger.debug('Line plot successfully generated')
        except Exception:
            logger.exception(f'Failed to generate line plot for {grand_id}')


if __name__ == '__main__':
    main()
