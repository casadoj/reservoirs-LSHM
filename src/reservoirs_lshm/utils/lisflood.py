from pathlib import Path
from datetime import datetime
from typing import Optional, Dict

import geopandas as gpd
import xarray as xr
import rioxarray as rio
from tqdm.auto import tqdm


def create_lisflood_nc(
    template: Path,
    points: gpd.GeoDataFrame,
    output: Optional[Path] = None,
    attrs: Optional[Dict] = None,
    crs: Optional[str] = None
) -> Optional[xr.DataArray]:
    """Creates a LISFLOOD-compatible NetCDF raster from a template and a set of points.

    This function loads a template NetCDF file, initializes it to zeros, and then
    assigns unique IDs from a GeoDataFrame to the nearest grid cells. It also
    handles adding metadata attributes and a Coordinate Reference System (CRS)
    before optionally saving the output to a new NetCDF file.

    Parameters
    ----------
    template : Path
        The file path to the template NetCDF file. This file must contain a 2D
        'template' data variable with 'lon' and 'lat' dimensions.
    points : geopandas.GeoDataFrame
        A GeoDataFrame where the index serves as the unique ID for each point
        and the 'geometry' column contains the point coordinates.
    output : Path, optional
        The file path to save the final NetCDF file. If not provided, the
        resulting xarray.DataArray is returned instead of being saved.
    attrs : dict, optional
        A dictionary of attributes to add to the output DataArray. The special
        key 'name' can be used to set the name of the DataArray.
    crs : str, optional
        A string representing the CRS (e.g., 'EPSG:4326') to write to the
        DataArray using the rioxarray accessor.

    Returns
    -------
    xarray.DataArray or None
        Returns the xarray.DataArray object if `output` is None, otherwise
        returns None.
    """
    # load template of the map
    ds = xr.open_dataset(template)['template']
    ds.close()

    # Special handling for 'name' from attrs
    if attrs and 'name' in attrs:
        ds.name = attrs['name']
    
    # remove all values
    ds[:,:] = 0
    
    # assign IDs
    for idx, point in tqdm(points.geometry.items(), total=len(points)):
        pixel = ds.sel({'lon': point.x, 'lat': point.y}, method='nearest')
        ds.loc[{'lon': pixel.lon.item(), 'lat': pixel.lat.item()}] = idx
    
    # modify attributes
    if attrs is not None:
        for key, value in attrs.items():
            if key != 'name':
                ds.attrs[key] = value
    
    # define coordinate reference system
    if crs is not None:
        ds = ds.rio.write_crs(crs)
    
    # export
    if output is not None:
        ds.to_netcdf(output)
        print(f'Created file {output}')
    else:
        return ds