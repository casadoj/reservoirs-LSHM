"""
Copyright 2019-2023 European Union
Licensed under the EUPL, Version 1.2 or as soon they will be approved by the European Commission  subsequent versions of the EUPL (the "Licence");
You may not use this work except in compliance with the Licence.
You may obtain a copy of the Licence at:
https://joinup.ec.europa.eu/sites/default/files/inline-files/EUPL%20v1_2%20EN(1).txt
Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the Licence for the specific language governing permissions and limitations under the Licence.

"""

import argparse
import os
from pathlib import Path
import sys
import time
from datetime import datetime
import xarray as xr
from typing import Dict, List, Union, Optional
from tqdm.auto import tqdm
import logging
import psutil
import traceback


def log_memory_usage():
    """Log memory usage to help detect issues related to segmentation faults."""
    process = psutil.Process()
    mem_info = process.memory_info()
    # Convert to MB
    logging.debug(f"Memory Usage: {mem_info.rss / (1024 * 1024):.2f} MB")


def check_coordinates(ds: Union[xr.Dataset, xr.DataArray]) -> Union[xr.Dataset, xr.DataArray]:
    """Makes sure that the geographical coordinates are named 'lat' and 'lon' and rounds the values to 9 decimals to avoid conflicts.
    """

    # check names of the coordinates
    ds = ds.rename(
        {k: v for k, v in {'x': 'lon', 'y': 'lat'}.items() if k in ds.sizes})

    # round it to avoid issues when comparing with other datasets
    ds['lon'] = ds['lon'].round(9)
    ds['lat'] = ds['lat'].round(9)

    return ds


def read_inputmaps(inputmaps: Union[str, Path]) -> xr.Dataset:
    """Reads NetCDF files from the input directory

    Parameters:
    -----------
    inputmaps: str or pathlib.Path
        directory that contains the input NetCDF files whose statistics will be computed. These files can be static (withouth time dimension) or dynamic (with time dimension)

    Returns:
    --------
    ds: xr.Dataset    
    """

    logging.info(f"Reading input maps from {inputmaps}")
    log_memory_usage()

    inputmaps = Path(inputmaps)
    if not inputmaps.is_dir():
        logging.error(f'{inputmaps} is missing or not a directory!')
        sys.exit(1)

    filepaths = list(inputmaps.glob('*.nc'))
    if not filepaths:
        logging.error(f'No NetCDF files found in {inputmaps}')
        sys.exit(2)

    logging.info(f'{len(filepaths)} NetCDF files found in {inputmaps}')

    try:
        # for dynamic maps
        ds = xr.open_mfdataset(
            filepaths,
            chunks='auto',
            # parallel=True, # ! the tool crashes when activated
            engine='netcdf4'
        )
        ds.close()
        # chunks is set to auto for general purpose processing
        # it could be optimized depending on input NetCDF
    except:
        # for static maps
        ds = xr.Dataset({file.stem.split('_')[0]: xr.open_dataset(
            file, engine='netcdf4')['Band1'] for file in filepaths})
        ds.close()

    if 'wgs_1984' in ds:
        ds = ds.drop_vars('wgs_1984')

    ds = check_coordinates(ds)

    logging.info("Successfully read input maps.")

    return ds


def read_masks(mask: Union[str, Path]) -> Dict[int, xr.DataArray]:
    """It loads the catchment masks in NetCDF formal from the input directory

    Parameters:
    -----------
    mask: str or pathlib.Path
        directory that contains the NetCDF files that define the catchment boundaries. These files can be the output of the `cutmaps` tool

    Returns:
    --------
    masks: dictionary of xr.DataArray
        keys represent the catchment ID and the values boolean maps of the catchment
    """

    logging.info(f"Reading masks from {mask}")
    log_memory_usage()

    # check masks
    mask = Path(mask)
    if not mask.is_dir():
        logging.error(f'{mask} is not a directory!')
        sys.exit(1)

    maskpaths = list(mask.glob('*.nc'))
    if not maskpaths:
        loggig.error(f'No NetCDF files found in "{mask}"')
        sys.exit(2)

    logging.info(f'{len(maskpaths)} mask NetCDF files found in "{mask}"')

    # load masks
    masks = {}
    for maskpath in maskpaths:
        ID = int(maskpath.stem)
        try:
            try:
                aoi = xr.open_dataset(maskpath, engine='netcdf4')['Band1']
                aoi.close()
            except:
                aoi = xr.open_dataarray(maskpath, engine='netcdf4')
                aoi.close()
            aoi = xr.where(aoi.notnull(), 1, aoi)
            aoi = check_coordinates(aoi)
            masks[ID] = aoi
        except Exception as e:
            logging.error(f'The mask {maskpath} could not be read: {e}')
            continue

    logging.info("Successfully read masks.")

    return masks


def read_pixarea(pixarea: Union[str, Path]) -> xr.DataArray:
    """It reads the LISFLOOD pixel area static map

    Parameters:
    -----------
    pixarea: string or Path
        a NetCDF file with pixel area used to compute weighted statistics. It is specifically meant for geographic projection systems where the area of a pixel varies with latitude

    Returns:
    --------
    weight: xr.DataArray
    """

    pixarea = Path(pixarea)
    if not pixarea.is_file():
        logging.error(f'{pixarea} is not a file!')
        sys.exit(1)

    try:
        weight = xr.open_dataset(pixarea, engine='netcdf4')['Band1']
        weight.close()
    except Exception as e:
        logging.error(f'The weighing map "{pixarea}" could not be loaded: {e}')
        sys.exit(2)

    # check coordinates
    weight = check_coordinates(weight)

    return weight


def catchment_statistics(maps: Union[xr.DataArray, xr.Dataset],
                         masks: Dict[int, xr.DataArray],
                         statistic: Union[str, List[str]],
                         weight: Optional[xr.DataArray] = None,
                         output: Optional[Union[str, Path]] = None,
                         overwrite: bool = False
                         ) -> Optional[xr.Dataset]:
    """
    Given a set of input maps and catchment masks, it computes catchment statistics.

    Parameters:
    -----------
    maps: xarray.DataArray or xarray.Dataset
        map or set of maps from which catchment statistics will be computed
    masks: dictionary of xr.DataArray
        a set of catchment masks. The tool `cutmaps` in this repository can be used to generate them
    statistic: string or list of strings
        statistics to be computed. Only some statistics are available: 'mean', 'sum', 'std', 'var', 'min', 'max', 'median', 'count'
    weight: optional or xr.DataArray
        map used to weight each pixel in "maps" before computing the statistics. It is meant to take into account the different pixel area in geographic projections
    output: optional, str or pathlib.Path
        directory where the resulting NetCDF files will be saved. If not provided, the results are put out as a xr.Dataset
    overwrite: boolean
        whether to overwrite or skip catchments whose output NetCDF file already exists. By default is False, so the catchment will be skipped

    Returns:
    --------
    A xr.Dataset of all catchment statistics or a NetCDF file for each catchment in the "masks" dictionary
    """

    logging.info("Starting catchment statistics computation...")
    log_memory_usage()

    start_time = time.perf_counter()

    if isinstance(maps, xr.DataArray):
        maps = xr.Dataset({maps.name: maps})

    # check statistic
    if isinstance(statistic, str):
        statistic = [statistic]
    possible_stats = ['mean', 'sum', 'std',
                      'var', 'min', 'max', 'median', 'count']
    assert all(stat in possible_stats for stat in statistic), "All values in 'statistic' should be one of these: {0}".format(
        ', '.join(possible_stats))
    stats_dict = {var: statistic for var in maps}

    # output directory
    if output is None:
        results = []
    else:
        output = Path(output)
        output.mkdir(parents=True, exist_ok=True)

    # define coordinates and variables of the resulting Dataset
    dims = dict(maps.sizes)
    dimnames = [dim.lower() for dim in dims]
    if 'lat' in dimnames and 'lon' in dimnames:
        x_dim, y_dim = 'lon', 'lat'
    else:
        x_dim, y_dim = 'x', 'y'
    del dims[x_dim]
    del dims[y_dim]
    coords = {dim: maps[dim] for dim in dims}
    variables = [f'{var}_{stat}' for var, stats in stats_dict.items()
                 for stat in stats]

    # compute statistics for each catchemnt
    for ID in tqdm(masks.keys(), desc='processing catchments'):
        # for ID in masks.keys():

        logging.info(f"Processing catchment ID {ID}")

        if output is not None:
            fileout = output / f'{ID:04}.nc'
            if fileout.exists() and not overwrite:
                logging.info(f"Skipping existing output file {fileout}")
                continue

        # create empty Dataset
        coords.update({'id': [ID]})
        maps_aoi = xr.Dataset(
            {var: xr.DataArray(coords=coords, dims=coords.keys()) for var in variables})

        # apply mask to the dataset
        aoi = masks[ID]
        masked_maps = maps.sel(
            {x_dim: aoi[x_dim], y_dim: aoi[y_dim]}).where(aoi == 1)
        masked_maps = masked_maps.compute()

        # apply weighting
        if weight is not None:
            masked_weight = weight.sel(
                {x_dim: aoi[x_dim], y_dim: aoi[y_dim]}).where(aoi == 1)
            weighted_maps = masked_maps.weighted(masked_weight.fillna(0))

        # compute statistics
        for var, stats in stats_dict.items():
            for stat in stats:
                if (stat in ['mean', 'sum', 'std', 'var']) and (weight is not None):
                    x = getattr(weighted_maps, stat)(dim=[x_dim, y_dim])[var]
                else:
                    x = getattr(masked_maps, stat)(dim=[x_dim, y_dim])[var]
                maps_aoi[f'{var}_{stat}'].loc[{'id': ID}] = x

        # save results
        if output is None:
            results.append(maps_aoi)
        else:
            maps_aoi.to_netcdf(fileout)

    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    logging.info(f"Processing completed in {elapsed_time:.2f} seconds")

    if output is None:
        results = xr.concat(results, dim='id')
        return results


def main(argv=sys.argv):
    prog = os.path.basename(argv[0])
    parser = argparse.ArgumentParser(
        description="""
        Utility to compute catchment statistics from (multiple) NetCDF files.
        The mask masp are NetCDF files with values in the area of interest and NaN elsewhere.
        The area map is optional and accounts for varying pixel area with latitude.
        """,
        prog=prog,
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Directory containing the input NetCDF files")
    parser.add_argument("-m", "--mask", required=True,
                        help="Directory containing the mask NetCDF files")
    parser.add_argument("-s", "--statistic", nargs='+', required=True,
                        help='List of statistics to be computed. Possible values: mean, sum, std, var, min, max, median, count')
    parser.add_argument("-o", "--output", required=True,
                        help="Directory where the output NetCDF files will be saved")
    parser.add_argument("-a", "--area", required=False, default=None,
                        help="NetCDF file of pixel area used to weigh the statistics")
    parser.add_argument("-w", "--overwrite", action="store_true",
                        default=False, help="Overwrite existing output files")

    args = parser.parse_args()

    # set up logging
    log_format = "%(asctime)s | %(levelname)s | %(message)s"
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    log_file = output_dir / f'catchstats_{datetime.now():%Y%m%d%H%M}.log'
    logging.basicConfig(
        level=logging.INFO,
        format=log_format,
        filename=log_file,
        filemode="w"
    )

    try:
        maps = read_inputmaps(args.input)
        masks = read_masks(args.mask)
        weight = read_pixarea(args.area) if args.area is not None else None
        catchment_statistics(maps, masks, args.statistic, weight=weight,
                             output=args.output, overwrite=args.overwrite)
    except Exception as e:
        logging.error(f'Fatal error: {e}')
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    sys.exit(main())
