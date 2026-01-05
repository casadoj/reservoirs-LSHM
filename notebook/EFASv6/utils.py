import numpy as np
import pandas as pd
import xarray as xr

def upstream_pixel(
    lat: float,
    lon: float,
    upArea: xr.DataArray
) -> tuple:
    """This function finds the upstream coordinates of a given point
    
    Parameteres:
    ------------
    lat: float
        latitude of the input point
    lon: float
        longitued of the input point
    upArea: xarray.DataArray
        map of upstream area
        
    Returns:
    --------
    lat: float
        latitude of the inmediate upstream pixel
    lon: float
        longitued of the inmediate upstream pixel
    """
    

    # Determine coordinate system
    lat_coord = 'lat' if 'lat' in upArea.coords else 'y'
    lon_coord = 'lon' if 'lon' in upArea.coords else 'x'
    
    try:

        # upstream area of the input coordinates
        area = upArea.sel({lat_coord: lat, lon_coord: lon}, method='nearest').item()
        
        # spatial resolution of the input map
        resolution = np.mean(np.diff(upArea[lon_coord].values))
        
        # Define window around the input pixel
        window = 1.5 * resolution
        upArea_window = upArea.sel({lat_coord: slice(lat + window, lat - window),
                                    lon_coord: slice(lon - window, lon + window)})
    
        # remove pixels with area smaller or equal than the input pixel
        mask = upArea_window.where(upArea_window < area, np.nan)

        # from the remaining pixels, find that with the largest upstream area
        up_pixel = mask.where(mask == mask.max(), drop=True)

        return (up_pixel[lat_coord].round(4).item(),
                up_pixel[lon_coord].round(4).item())
                
    except Exception as e:
        print(f"Failed to find upstream pixel: {e}")
        return None, None

def downstream_pixel(
    lat: float,
    lon: float,
    upArea: xr.DataArray
) -> tuple:
    """
    This function finds the downstream pixel coordinates of a given point.
    
    Parameters:
    -----------
    lat: float
        Latitude of the input point.
    lon: float
        Longitude of the input point.
    upArea: xarray.DataArray
        Map of upstream area.
        
    Returns:
    --------
    tuple:
        Latitude and longitude of the immediate downstream pixel,
        or (None, None) if not found.
    """
    
    # Determine coordinate system
    lat_coord = 'lat' if 'lat' in upArea.coords else 'y'
    lon_coord = 'lon' if 'lon' in upArea.coords else 'x'
    
    try:
        # upstream area of the input coordinates
        area = upArea.sel({lat_coord: lat, lon_coord: lon}, method='nearest').values
        
        # spatial resolution of the input map
        resolution = np.mean(np.diff(upArea[lon_coord].values))
        
        # Define window around the input pixel
        window = 1.5 * resolution #np.array([-1.5 * resolution, 1.5 * resolution])
        upArea_window = upArea.sel({lat_coord: slice(lat + window, lat - window),
                                    lon_coord: slice(lon - window, lon + window)})
        
        # Mask out pixels with area less than or equal to the input pixel
        mask = upArea_window.where(upArea_window > area, drop=True)
        
        # Find the pixel with the smallest upstream area
        down_pixel = mask.where(mask == mask.max(), drop=True)
        
        return (down_pixel[lat_coord].round(4).item(),
                down_pixel[lon_coord].round(4).item())
                
    except Exception as e:
        print(f"Failed to find downstream pixel: {e}")
        return None, None

def outlet_width(
    chanbw: xr.DataArray,
    uparea: xr.DataArray,
    x: float,
    y: float,
    n_points: int = 5):
    """It searches for the width of the outlet of a water body. Given the coordinates of the lake/reservoir in the LISFLOOD grid, it finds the 'n_points' upstream and downstream, and it selects the width of the pixel with a largest change in the sequence.
    
    Parameters:
    -----------
    chanbw: xarray.DataArray
        static map of channel bankfull width
    uparea: xarray.DataArray
        static map of upstream area
    x: float
        coordinate X of the estimated outlet of the water body
    y: float
        coordinate Y of the estimated outlet of the water body
    n_points: int
        number of pixels to search. This number will be applied both upstream and downstream
        
    Returns:
    --------
    width: float
        width of the water body outlet
    """
    
    xs = np.array([x])
    ys = np.array([y])
    
    n_points = 2
    # coordinates of the two pixels upstream and downstream
    i = 0
    while i < n_points:
        try:
            y_up, x_up = upstream_pixel(ys[0], xs[0], uparea)
            xs = np.hstack([x_up, xs])
            ys = np.hstack([y_up, ys])
        except:
            break
        i += 1
        
    # coordinates of the two pixels downstream
    i = 0
    while i < n_points:
        try:
            y_down, x_down = downstream_pixel(ys[-1], xs[-1], uparea)
            if (y_down is not None) and (x_down is not None):
                xs = np.hstack([xs, x_down])
                ys = np.hstack([ys, y_down])
        except:
            break
        i += 1
        
    # extract channel width
    width = np.array([chanbw.sel(x=x, y=y, method='nearest').item() for x, y in zip(xs, ys)])
    
    return np.nanmin(width)