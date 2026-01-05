#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import logging
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm
from datetime import datetime

from . import Config, read_attributes, read_timeseries
from .models.starfit.storage import fit_storage, create_storage_harmonic
from .models.starfit.release import fit_release, create_release_harmonic, create_release_linear
from .models.starfit.functions import plot_nor, plot_release
from .utils.logging import setup_logger


def main():

    # parse arguments
    parser = argparse.ArgumentParser(
        description="""
        Fit the storage and release rules for the Starfit reservoir routine.
        The fitted models are saved as Pickle files and plotted against the
        observed data used for fitting.
        """
        )
    parser.add_argument('-c', '--config-file', type=str, required=True, help='Path to the configuration file')
    parser.add_argument('-o', '--overwrite', action='store_true', default=False, help='Overwrite existing model.')
    args = parser.parse_args()

    # set up logger
    logger = setup_logger(
        name=__name__,
        log_level=logging.INFO,
        log_file=f'{datetime.now():%Y%m%d%H%M%S}_fit_starfit.log'
    )

    # read configuration file
    cfg = Config(args.config_file)
    PATH_STORAGE = cfg.PATH_DEF.parent / 'NOR'
    PATH_RELEASE = cfg.PATH_DEF.parent / 'release'
    PATH_STORAGE.mkdir(parents=True, exist_ok=True)
    PATH_RELEASE.mkdir(parents=True, exist_ok=True)
    logger.info(f'Storage models will be saved in: {PATH_STORAGE}')
    logger.info(f'Release models will be saved in: {PATH_RELEASE}')

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
        for grand_id, obs in timeseries.items():
            # convert units
            obs['s'] = obs.storage * 1e-6 # MCM
            obs[['i', 'r']] = obs[['inflow', 'outflow']] * 1e-6 * 86400 # MCM/day
            # update reservoir capacity, if maximum observation exceeds GRanD
            attributes.loc[grand_id, 'CAP_MCM'] = max(attributes.loc[grand_id, 'CAP_MCM'], obs.s.max())
        logger.info(f'{len(timeseries)} reservoirs with timeseries')
    except IOError:
        logger.exception(f'Failed to read time series')
        raise

    # === Fit Storage Models ===
    for grand_id, obs in tqdm(timeseries.items(), desc="Fitting storage models"):
        
        out_file = PATH_STORAGE / f'{grand_id}.pkl'
        if out_file.exists() and not args.overwrite:
            logger.info(f'Storage model already exists for {grand_id}, skipping (use --overwrite to force)')
            continue

        logger.info(f'Fitting storage model for {grand_id}')

        # fit storage model
        try:
            for years, n_points in zip([8, 6, 4], [3, 2, 2]):
                model_storage = fit_storage(
                    grand_id,
                    storage_daily=obs.s,
                    attributes=attributes.loc[grand_id],
                    min_days=years * 365,
                    n_points=n_points,
                )
                if not model_storage['weekly storage'].empty:
                    break
    
            if model_storage['weekly storage'].empty:
                logger.warning(f'Could not fit storage model for {grand_id}')
                continue
                
        except Exception:
            logger.exception(f'Failed to fit storage model for {grand_id}')
            continue

        # export fitted parameters
        try:
            pars_to_export = ['capacity (MCM)', 'NOR upper bound', 'NOR lower bound']
            pars = {key: value for key, value in model_storage.items() if key in pars_to_export}
            with open(out_file, 'wb') as file:
                pickle.dump(pars, file)
            logger.debug(f'Saved storage model for {grand_id}')
        except Exception as e:
            logger.exception(f'Failed to save the storage model for {grand_id}')

        # define normal operating range (NOR)
        try:
            NORup = create_storage_harmonic(model_storage['NOR upper bound'], name='flood').set_index('epiweek')
            NORdown = create_storage_harmonic(model_storage['NOR lower bound'], name='conservation').set_index('epiweek')
            NOR = pd.concat((NORup, NORdown), axis=1)
        except Exception as e:
            logger.exception(f'Failed to create the normal operating rules (NOR) for {grand_id}')
    
        # weekly time series of standardised storage combined with NOR
        weekly_storage = model_storage['weekly storage']
    
        # plot model
        try:
            plot_nor(
                weekly_storage,
                NOR,
                n_points=n_points,
                title='{0} - {1}'.format(grand_id, attributes.loc[grand_id, 'DAM_NAME']),
                save=PATH_STORAGE / f'{grand_id}.jpg'
            )
        except Exception:
            logger.exception(f'Failed to plot the storage model of reservoir {grand_id}')

    # === Fit Release Models ===
    grand_ids = [int(file.stem) for file in PATH_STORAGE.glob('*.pkl')]
    for grand_id, obs in tqdm(timeseries.items(), desc="Fitting storage models"):

        if grand_id not in grand_ids:
            logger.info(f"Skipping {grand_id} as it doesn't have a storage model")
            
        out_file = PATH_RELEASE / f'{grand_id}.pkl'
        if out_file.exists() and not args.overwrite:
            logger.info(f'Release model already exists for {grand_id}, skipping (use --overwrite to force)')
            continue

        logger.info(f'Fitting release model for {grand_id}')

        # fit release model
        try:
            for years in [5, 4]:
                model_release = fit_release(
                    grand_id,
                    daily_ops=obs[['s', 'i', 'r']],
                    attributes=attributes.loc[grand_id],
                    NOR_path=PATH_STORAGE,
                    cutoff_year=None,
                    min_weeks=52 * years
                )
                if pd.notna(model_release['mean inflow (MCM/wk)']):
                    break

            if not model_release or all(np.isnan(model_release['harmonic parameters'])):
                logger.warning(f'Could not fit release model for {grand_id}')
                continue

        except Exception:
            logger.exception(f'Failed to fit release model for {grand_id}')

        # export fitted parameters
        try:
            pars_to_export = ['mean inflow (MCM/wk)', 'harmonic parameters', 'residual parameters', 'constraints']
            pars = {key: value for key, value in model_release.items() if key in pars_to_export}
            with open(out_file, 'wb') as file:
                pickle.dump(pars, file)
        except Exception:
            logger.exception(f'Failed to save the release model for {grand_id}')

        # extract info from the fitted release: average inflow, harmonic release (standardised) and release contraints
        try:
            avg_inflow = model_release['mean inflow (MCM/wk)']
            release_harmonic = create_release_harmonic(model_release['harmonic parameters']).set_index('epiweek').squeeze()
            release_linear = create_release_linear(model_release['residual parameters'])
            Qmin, Qmax = model_release['constraints']
            weekly_release = model_release['weekly release'].set_index('epiweek')
        except Exception as e:
            logger.exception(f'Failed to create the release rules for {grand_id}')
    
        # plot model
        try:
            title = '{0} - {1}'.format(grand_id, attributes.loc[grand_id, 'DAM_NAME'])
            plot_release(
                weekly_release.r, 
                avg_inflow, 
                release_harmonic, 
                release_linear, 
                Qmin, 
                Qmax, 
                title=title,
                save=PATH_RELEASE / f'{grand_id}.jpg'
            )
        except Exception:
            logger.exception(f'Failed to plot the release model for {grand_id}')


if __name__ == "__main__":
    main()
