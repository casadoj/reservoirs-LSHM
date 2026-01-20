import numpy as np
import pandas as pd
from scipy.interpolate import interp1d, PchipInterpolator
from scipy.stats import gaussian_kde
from typing import Literal, Union, Optional
import matplotlib.pyplot as plt


def bin_data(
    elevation: pd.Series, 
    target: Union[pd.Series, pd.DataFrame], 
    agg: Literal['median', 'mean'] = 'median',
    bin_size: float = 0.5,
    ) -> pd.Series:
    """
    Bins reservoir elevation and corresponding storage data into regular elevation intervals
    and computes the mean storage for each bin.

    Parameters
    ----------
    elevation : pd.Series
        Series of elevation values (in meters), typically from time series data.
    target: Union[pd.Series, pd.DataFrame]
        Series of storage, area or other variable corresponding to the elevation series.
    agg: Literal['median', 'mean']
        Statistic used to bin the input data
    bin_size : float, optional
        The elevation bin size (in meters) to aggregate the data, default is 0.1 m.

    Returns
    -------
    pd.Series
        Series with binned elevation values as the index and the mean storage for each bin.
        The index represents the center of each elevation bin.
    """

    if isinstance(target, pd.Series):
        target_df = pd.DataFrame(target)
    else:
        target_df = target.copy()
    df = pd.concat([elevation.rename('elevation'), target_df], axis=1).dropna(axis=1, how='any')
    df.sort_values('elevation', inplace=True)
    df.reset_index(drop=True, inplace=True)

    # Define bins: from min to max elevation, spaced every bin_size
    min_elev = np.ceil(df.elevation.min() / bin_size) * bin_size
    max_elev = np.floor(df.elevation.max() / bin_size) * bin_size
    bins = np.round(np.arange(min_elev, max_elev + .01, bin_size), 3)
    #bins = np.append(np.append(df.elevation.min(), bins), df.elevation.max())
        
    # bin the elevation values
    df['elev_bin'] = pd.cut(df.elevation, bins, include_lowest=False)

    # group by bin and compute mean storage (and optionally elevation)
    agg_dict = {col: agg for col in target_df.columns}
    binned = df.groupby('elev_bin', observed=False).agg(agg_dict)

    # replace bin labels with bin centers
    binned.index = np.mean([bins[:-1], bins[1:]], axis=0)
    binned.index.name = 'elevation'

    # remove bins with no data
    binned.dropna(how='any', inplace=True)
        
    if any(binned.diff().min() < 0):
        print('WARNING. The binned data is not monotonically increasing')

    return binned.squeeze()


def remove_outliers_kde(
        df: pd.DataFrame, 
        elevation_col: str = 'elevation', 
        storage_col: str = 'storage', 
        threshold_density: float = 0.005,
        inplace: bool = False
    ) -> pd.DataFrame:
    """
    Removes outliers from a 2D scatter plot (e.g., elevation vs. storage) 
    based on Kernel Density Estimation (KDE) thresholding.

    This method identifies sparse points (outliers) in the 2D distribution 
    by calculating the probability density at each point.

    Args:
        df (pd.DataFrame): The input DataFrame.
        elevation_col (str): The name of the column containing the first variable (x-axis), 
                             default is 'elevation'.
        storage_col (str): The name of the column containing the second variable (y-axis), 
                           default is 'storage'.
        threshold_density (float): The minimum density value (KDE output) for a point to be 
                                   considered an 'inlier'. Points with density below this 
                                   value are removed. This value often needs manual tuning.
        inplace (bool): If True, the original DataFrame 'df' is modified by dropping the 
                        outlier rows. If False (default), a new DataFrame with only inliers 
                        is returned.

    Returns:
        pd.DataFrame or None: If 'inplace' is False, a new DataFrame containing only the 
                              inlier data points is returned. If 'inplace' is True, the 
                              original DataFrame is modified and None is returned.
    """
    # 1. Prepare Data
    notnan_mask = df[[elevation_col, storage_col]].notna().all(axis=1)
    array = df[notnan_mask][[elevation_col, storage_col]].T.values
    
    # 2. Perform 2D Kernel Density Estimation (KDE)
    kde = gaussian_kde(array)

    # 3. Evaluate the KDE for every point
    density = kde(array)
    inlier_mask = density >= threshold_density

    # 4. Filter the DataFrame using the mask
    inlier_df = df[notnan_mask][inlier_mask].copy()
    
    # 5. Report and Return
    num_outliers = array.shape[1] - len(inlier_df)
    print(f"Points retained (inliers): {len(inlier_df)}")
    print(f"Points removed (outliers): {num_outliers}")

    if inplace:
        df.drop(df.index.difference(inlier_df.index), inplace=True)
        return None
    else:
        return inlier_df

    
def fit_reservoir_curve(
    x_binned: pd.Series,
    y_binned: pd.Series, 
    method: Literal['poly1d', 'interp1d', 'pchip'] = 'pchip',
    degree: int = 2
):
    """
    Fits a smooth curve to a binned elevation-storage series using a selected method.

    This function models the relationship between reservoir elevation and storage 
    by fitting a smooth curve to binned data (typically pre-processed using fixed 
    elevation intervals). It supports polynomial fitting, linear interpolation, 
    and shape-preserving cubic Hermite interpolation (PCHIP).

    Parameters
    ----------
    x_binned : pd.Series
        A pandas Series representing the explanatory variable in the reservoir curve (e.g. elevation)
    y_binned: pd.Series
        A pandas Series representing the variable to be inferred with the reservoir curve (e.g. storage)
    method : {'poly1d', 'interp1d', 'pchip'}, optional
        The fitting method to use:
        - 'poly1d' fits a polynomial of specified degree (e.g., quadratic).
        - 'interp1d' performs linear interpolation.
        - 'pchip' uses shape-preserving cubic Hermite interpolation (default).
    degree : int, optional
        Degree of the polynomial if `method='poly1d'`. Ignored for other methods. 
        Default is 2.

    Returns
    -------
    callable
        A function that takes elevation values as input and returns estimated 
        storage values. The return type depends on the method:
        - `np.poly1d` for polynomial fitting,
        - `scipy.interpolate.interp1d` for linear interpolation,
        - `scipy.interpolate.PchipInterpolator` for PCHIP.
    
    Raises
    ------
    ValueError
        If an unsupported fitting method is specified.
    """

    if method.lower() == 'poly1d':
        coefficients = np.polyfit(x_binned, y_binned, degree)
        reservoir_curve = np.poly1d(coefficients)
    elif method.lower() == 'interp1d':
        reservoir_curve = interp1d(
            x=x_binned,
            y=y_binned,
            kind='linear',
            fill_value='extrapolate',
            assume_sorted=True
            )
    elif method.lower() == 'pchip':
        reservoir_curve = PchipInterpolator(
            x=x_binned,
            y=y_binned
        )
    else:
        raise ValueError(f'"method" must be either "interp1d" or "pchip": {method} was provided')

    return reservoir_curve


def storage_from_elevation(
    reservoir_curve: callable,
    elevation: Union[pd.Series, np.ndarray]
) -> Union[pd.Series, np.ndarray]:
    """
    Produces a time series of reservoir storage given the reservoir curve and an elevation time series.

    Parameters:
    -----------
    reservoir_curve: callable
        A NumPy polynomial object representing a fitted reservoir curve (storage vs elevation)
    elevation: pandas.Series or numpy.ndarray
        Reservoir elevation data

    Returns:
    --------
    storage: pandas.Series or numpy.ndarray
        Estimated reservoir storage data.
    """

    # estimate storage
    storage = reservoir_curve(elevation)
    if isinstance(elevation, pd.Series):
        storage = pd.Series(
            data=storage,
            index=elevation.index,
            name='storage'
            )

    return storage


def elevation_from_storage(
    reservoir_curve: np.poly1d,
    storage: pd.Series
) -> pd.Series:
    """
    Produces a time series of reservoir elevation given the reservoir curve and a storage time series.

    Parameters:
    -----------
    reservoir_curve: numpy.poly1d
        A NumPy polynomial object representing a fitted reservoir curve (storage vs elevation)
    storage: pandas.Series
        A pandas Series containing corresponding reservoir storage data.

    Returns:
    --------
    elevation: pandas.Series
        A pandas Series containing elevation data.
    """

    # coefficients of the polynomial
    a, b, c = reservoir_curve.coefficients

    # estimate elevation
    elevation = pd.Series(
        data=(-b + np.sqrt(b**2 - 4 * a * (c - storage))) / (2 * a),
        index=storage.index,
        name='elevation'
    )

    return elevation

def area_from_elevation(
    reservoir_curve: np.poly1d,
    elevation: pd.Series
) -> pd.Series:
    """
    Produces a time series of reservoir area given the reservoir curve and an elevation time series.

    The derivatie of the reservoir curve (storage-elevation) is the area-elevation curve:

            V = f(Z)

            A = dV / dZ = f'(Z)

    Parameters:
    -----------
    reservoir_curve: numpy.poly1d
        A NumPy polynomial object representing a fitted reservoir curve (storage vs elevation)
    elevation: pandas.Series
        A pandas Series containing elevation data.

    Returns:
    --------
    area: pandas.Series
        A pandas Series containing corresponding reservoir area data.
    """

    # estimate area
    try:
        area = pd.Series(
            data=reservoir_curve.deriv()(elevation),
            index=elevation.index,
            name='area'
            )
    except:
        area = pd.Series(
            data=reservoir_curve.derivative()(elevation),
            index=elevation.index,
            name='area'
        )
        
    return area


def elevation_sequence(
    z_min: float, 
    z_max: float, 
    method: Literal['linear', 'cosine', 'arctanh'] = 'linear', 
    step: int = 1, 
    N: int = 25, 
    alpha: float = .95
) -> np.ndarray:
    """
    Generates a sequence of elevation values within a specified range 
    using various spacing methods to control point density.

    Args:
        z_min (float): The minimum elevation value.
        z_max (float): The maximum elevation valu (end of the sequence).
        method (Literal['linear', 'cosine', 'arctanh'], optional): The method
            used for spacing the points:
            - 'linear': Uniform spacing defined by 'step'.
            - 'cosine': Clusters points at the **extremes** (z_min and z_max).
            - 'arctanh': Clusters points in the **middle** of the range.
            Defaults to 'linear'.
        step (int, optional): The step size used when `method` is 'linear'. 
            Ignored for 'cosine' and 'arctanh' methods. Defaults to 1.
        N (int, optional): The total number of points in the generated sequence
            when `method` is 'cosine' or 'arctanh'. Ignored for 'linear'. 
            Defaults to 25.
        alpha (float, optional): The clustering factor for the 'arctanh' method. 
            Must be less than 1 (e.g., 0.9 to 0.999). A higher value creates 
            tighter clustering (smaller steps) in the center. Ignored for 
            'linear' and 'cosine'. Defaults to 0.95.

    Returns:
        np.ndarray: A strictly increasing array of elevation values.

    Raises:
        ValueError: If an unrecognized string is passed to the 'method' parameter.
    """
    dam_hgt_m = z_max - z_min
    if method == 'linear':
        z_values = np.arange(z_min, z_max + .1, step)
        if z_max not in z_values:
            z_valuess = np.append(z_values, z_max)
    elif method == 'cosine':
        i = np.linspace(0, 1, N)
        clustered = 0.5 * (1 - np.cos(np.pi * i))
        z_values = z_min + dam_hgt_m * clustered
    elif method == 'arctanh':
        i = np.linspace(-1, 1, N)
        clustered = np.arctanh(i * alpha)
        z_values = z_min + dam_hgt_m * (clustered - clustered.min()) / (clustered.max() - clustered.min())
    else:
        raise ValueError(f"Method '{method}' not recognized. Must be 'linear', 'cosine', or 'arctanh'.")
    return z_values


def estimate_area_curve(
        lookup_table: pd.DataFrame, 
        elevation_col: str = 'elevation', 
        storage_col: str = 'storage'
    ) -> pd.Series:
    """
    Estimate area curve from elevation and storage data.
    
    Parameters:
    - lookup_table: DataFrame with elevation and storage columns.
    - elevation_col: Name of the elevation column.
    - storage_col: Name of the storage column.
    
    Returns:
    - Series of the area associated to the entries in the lookup table.
    """
    area = lookup_table[storage_col].diff() / lookup_table[elevation_col].diff()
    if lookup_table[storage_col].iloc[0] == 0:
        for i, idx in enumerate(lookup_table.index):
            if i == 0:
                area.loc[idx] = 0
            else:
                area.loc[idx] = 2 * area.iloc[i] - area.iloc[i - 1]

    return area


def plot_reservoir_curves(
    reservoir_curve: Optional[pd.DataFrame] = None,
    attrs: Optional[pd.Series] = None,
    obs: Optional[pd.DataFrame] = None,
    **kwargs
    ):
    """
    Generates a 2x2 matrix of scatter plots showing the relationships between
    reservoir elevation, area, and storage (volume), which are collectively 
    known as the reservoir's characteristic curves.

    The function plots the three essential relationships:
    1. Elevation vs. Storage (Volume)
    2. Elevation vs. Area
    3. Area vs. Storage (Volume)

    The fourth subplot (Area vs. Area) is left blank. Reference lines for 
    key values (e.g., minimum/maximum elevation, maximum area/storage) are 
    drawn based on external variables (`elev_masl`, `dam_hgt_m`, `area_skm`, 
    `cap_mcm`) that must be defined in the global or enclosing scope.

    Args:
        reservoir_curve (pd.DataFrame, optional): A DataFrame containing the reservoir 
            characteristic curve data. It must contain the following columns:
            - 'elevation' (float): Reservoir elevation (masl).
            - 'area' (float): Reservoir surface area (km2).
            - 'volume' or 'storage' (float): Reservoir storage volume (hm3).
        attrs (pd.Series, optional): An optional Series containing reservoir
            attributes from GRanD or GDW such as 'DAM_HGT_M', 'ELEV_MASL', 
            'AREA_SKM' or 'CAP_MCM' to draw reference lines on the plots. 
            Defaults to None.
        obs (pd.DataFrame, optional): An optional DataFrame containing observed
            reservoir data to overlay on the plots. It should have the same
            columns as `reservoir_curve` for the variables being plotted.
            Defaults to None.
            
        **kwargs: Optional keyword arguments to customize the plot:
            - figsize (tuple, optional): Size of the figure (width, height). 
              Defaults to (10, 10).

    Returns:
        tuple[plt.Figure, np.ndarray]: A tuple containing:
            - fig (matplotlib.figure.Figure): The main Matplotlib figure object.
            - axes (np.ndarray): A 2x2 array of Matplotlib axes objects.

    Notes:
        This function relies on external variables for reference lines that 
        must be accessible in the function's scope, including:
        - `elev_masl`: Maximum elevation (masl).
        - `dam_hgt_m`: Dam height (m) or relative minimum elevation.
        - `area_skm`: Maximum area (km2).
        - `cap_mcm`: Maximum storage capacity (hm3).
    """
    alpha = kwargs.get('alpha', 0.3)
    cmap = kwargs.get('cmap', 'coolwarm')
    figsize = kwargs.get('figsize', (10, 10))
    size = kwargs.get('size', 8)
    

    fig, axes = plt.subplots(ncols=2, nrows=2, figsize=figsize, sharex='col', sharey='row')

    if attrs is not None:
        dam_hgt_m, elev_masl, cap_mcm, area_skm = attrs.loc[['DAM_HGT_M', 'ELEV_MASL', 'CAP_MCM', 'AREA_SKM']]
    var_props = {
        'elevation': {
            'label': 'elevation (masl)',
            'ref': [elev_masl - dam_hgt_m, elev_masl] if attrs is not None else []
        },
        'area': {
            'label': 'area (km2)',
            'ref': [0, area_skm] if attrs is not None else [0]
        },
        'storage': {
            'label': 'volume (hm3)',
            'ref': [0, cap_mcm] if attrs is not None else [0]
        }
    }

    aux_props = dict(ls='--', lw=.5, c='k', zorder=0)
    curve_props = dict(lw=1, c='k', zorder=2)

    for j, var_x in enumerate(['elevation', 'area']):
        for i, var_y in enumerate(['storage', 'area']):
            
            ax = axes[i,j]
            if i == 1 & j == 1:
                ax.axis('off')
                continue
                
            if reservoir_curve is not None:
                if var_x == 'elevation':
                    if var_y == 'storage':
                        ax.plot(reservoir_curve.elevation, reservoir_curve.volume, **curve_props)
                    if var_y == 'area':
                        ax.plot(reservoir_curve.elevation, reservoir_curve.area, **curve_props)
                elif var_x == 'area':
                    if var_y == 'storage':
                        ax.plot(reservoir_curve.area, reservoir_curve.volume, **curve_props)

            # scatter plot of observed data
            if obs is not None and all(col in obs.columns for col in [var_x, var_y]):
                ax.scatter(
                    obs[var_x], 
                    obs[var_y], 
                    c=obs.index, 
                    cmap=cmap,
                    s=size,
                    alpha=alpha,
                    zorder=1,
                    label='observations'
                )
                
            for x in var_props[var_x]['ref']:
                ax.axvline(x, **aux_props)
            for y in var_props[var_y]['ref']:
                ax.axhline(y, **aux_props)

            if (i == 1) | (j == 1):
                ax.set_xlabel(var_props[var_x]['label'])
            if j == 0:
                ax.set_ylabel(var_props[var_y]['label'])

    return fig, axes