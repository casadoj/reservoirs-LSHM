import numpy as np
import pandas as pd
from scipy.interpolate import interp1d, PchipInterpolator
from scipy.stats import gaussian_kde
from typing import Literal, Union, Optional
import matplotlib.pyplot as plt
import pickle
from pathlib import Path


def bin_data(
    elevation: pd.Series, 
    target: Union[pd.Series, pd.DataFrame], 
    agg: Literal['median', 'mean', 'closest'] = 'median',
    bin_size: float = 0.5,
    ) -> Union[pd.Series, pd.DataFrame]:
    """
    Bins reservoir elevation and corresponding storage data into regular elevation intervals.

    Parameters
    ----------
    elevation : pd.Series
        Series of elevation values (in meters), typically from time series data.
    target: Union[pd.Series, pd.DataFrame]
        Series of storage, area or other variable corresponding to the elevation series.
    agg: Literal['median', 'mean'. 'closest']
        Statistic used to bin the input data. If 'closest', the closest observation to each bin center is used.
        Default is 'median'.
    bin_size : float, optional
        The elevation bin size (in meters) to aggregate the data, default is 0.5 m.

    Returns
    -------
    Union[pd.Series, pd.DataFrame]
        Series with binned elevation values as the index and the mean storage for each bin.
        The index represents the center of each elevation bin.
    """

    if isinstance(target, pd.Series):
        target_df = pd.DataFrame(target)
    else:
        target_df = target.copy()
    df = pd.concat([elevation.rename('elevation'), target_df], axis=1).dropna(axis=1, how='any')
    df.drop_duplicates(inplace=True)
    df.sort_values('elevation', inplace=True)
    df.reset_index(drop=True, inplace=True)

    # Define bins: from min to max elevation, spaced every bin_size
    min_elev = np.ceil(df.elevation.min() / bin_size) * bin_size
    max_elev = np.floor(df.elevation.max() / bin_size) * bin_size
    bins = np.round(np.arange(min_elev, max_elev + bin_size / 10, bin_size), 3)
        
    if agg == 'closest':
        # add minimum and maximum observed elevation
        #bins = np.append(np.append(df.elevation.min(), bins), df.elevation.max())
        
        # keep only the closest observation to each bin
        diff_matrix = np.abs(bins[:, np.newaxis] - df.elevation.values)
        closest_indices = np.argmin(diff_matrix, axis=1)
        binned = df.iloc[closest_indices].copy()
        #binned.set_index('elevation', inplace=True, drop=True)
        binned.reset_index(drop=True, inplace=True)
        
    elif agg in ['mean', 'median']:
        # bin the elevation values
        df['elev_bin'] = pd.cut(df.elevation, bins, include_lowest=False)

        # group by bin and compute mean storage (and optionally elevation)
        agg_dict = {col: agg for col in target_df.columns}
        binned = df.groupby('elev_bin', observed=False).agg(agg_dict)

        # replace bin labels with bin centers
        binned.index = np.mean([bins[:-1], bins[1:]], axis=0)
        binned.index.name = 'elevation'
        binned.reset_index(inplace=True)

        # remove bins with no data
        binned.dropna(how='any', inplace=True)

    else:
        raise ValueError(f'"agg" must be either "median", "mean" or "closest": {agg} was provided')
        
    if any(binned.diff().min() < 0):
        print('WARNING. The binned data is not monotonically increasing')

    return binned


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


class ReservoirCurve(pd.DataFrame):
    """
    Enhanced pandas.DataFrame subclass for Elevation-Area-Storage (EAS) curve analysis.

    The ReservoirCurve object stores the lookup table data and provides methods
    for fitting and inferring values between elevation, storage, and area time series.
    It automatically enforces monotonicity and respects defined physical limits
    (z_min, z_max, v_max, etc.) to prevent non-physical extrapolation.

    Attributes:
        z_min (float): Minimum physical elevation (e.g., dam invert).
        z_max (float): Maximum physical elevation (e.g., dam crest).
        v_min (float): Minimum observed storage (always >= 0).
        v_max (float): Maximum physical storage capacity.
        curve_zv (callable): Fitted interpolator for Elevation -> Storage (set by .fit()).
        curve_vz (callable): Fitted interpolator for Storage -> Elevation (set by .fit()).
    """

    _metadata = ['z_min', 'z_max', 'v_min', 'v_max', 'curve_zv', 'curve_vz']
    
    def __init__(
            self, 
            lookup_table: pd.DataFrame,
            z_min: Optional[float] = None,
            z_max: Optional[float] = None,
            v_min: Optional[float] = None,
            v_max: Optional[float] = None,
            *args, 
            **kwargs
            ):
        """
        Initializes the ReservoirCurve object.

        Parameters
        ----------
        lookup_table : pd.DataFrame
            A DataFrame containing the EAS curve data. Must include 'elevation' (z)
            and 'storage' (v). Data must be monotonically
            increasing with elevation.
        *args, **kwargs :
            Arguments passed to the pandas.DataFrame constructor.

        Raises
        ------
        ValueError
            If 'elevation' or 'storage' columns are missing, or if storage/area
            are not monotonically increasing with elevation, or if observed data
            exceeds user-defined limits.
        """

        super().__init__(lookup_table, *args, **kwargs)
        self.sort_values('elevation', inplace=True)
        self.reset_index(inplace=True, drop=True)

        # check monotonicity in the storage values
        if not self['storage'].is_monotonic_increasing:
            raise ValueError("The 'storage' column must be monotonically increasing with 'elevation'. Check data quality.")
        
        # estimate area
        if 'area' not in lookup_table.columns:
            self['area'], curve_za = self._lookup_area()
        # check monotonicity in the area values
        if not self['area'].is_monotonic_increasing:
            raise ValueError("The 'area' column must be monotonically increasing with 'elevation'. Check data quality.")

        # define curve limits
        self.z_min = z_min if z_min is not None else self['elevation'].min()
        self.z_max = z_max if z_max is not None else self['elevation'].max()
        self.v_min = v_min if v_min is not None else self['storage'].min()
        self.v_max = v_max if v_max is not None else self['storage'].max()
        self.a_min = curve_za(self.z_min)
        self.a_max = curve_za(self.z_max)

        # initialize empty curves
        self.curve_zv = None
        self.curve_vz = None
        self.curve_za = None
        self.curve_az = None
        self.curve_av = None
        self.curve_va = None

    def _check_range(self, data: Union[pd.Series, np.ndarray], variable: Literal['elevation', 'storage']) -> Union[pd.Series, np.ndarray]:
        """Converts into NaN values outside the reservoir curve range to avoid extrapolation problems

        Parameters:
        -----------
        data: pandas.Series or numpy.ndarray
            Values to be checked
        variable: string
            Defines the variable of "data"

        Returns:
        --------
        np.ndarray
            The input data with out-of-range values set to NaN.
        """
        array = np.array(data)
        
        if variable == 'elevation':
            min_value, max_value = self.z_min, self.z_max
        elif variable == 'storage':
            min_value, max_value = self.v_min, self.v_max
        elif variable == 'area':
            min_value, max_value = self.a_min, self.a_max
        else:
            raise ValueError(f'"variable" must be either "elevation", "storage" or "area": {variable} was provided')
        
        mask = (array < min_value) | (array > max_value)
        if mask.sum() > 0:
            array[mask] = np.nan
            print(f'WARNING. {mask.sum()} {variable} values were removed because they were outside of the range [{min_value:.3f},{max_value:.3f}]')

        return array
    
    def _convert(
        self, 
        input: Union[pd.Series, np.ndarray], 
        curve_attr: str,
        input_var: Literal['elevation', 'storage', 'area'], 
        output_var: Literal['elevation', 'storage', 'area']
    ) -> Union[pd.Series, np.ndarray]:
        """Generic method to convert a time series using a fitted curve."""
        
        # check input data is within the allowed range
        input = input.copy()
        self._check_range(input, variable=input_var)

        # get the fitted curve
        curve = getattr(self, curve_attr)
        if curve is None:
             raise ValueError(f"The curve '{curve_attr}' has not been fitted. Call .fit() first.")

        # convert data and check range
        output = curve(input)
        self._check_range(output, variable=output_var)
        if isinstance(input, pd.Series):
            output = pd.Series(data=output, index=input.index, name=output_var)
        
        return output
    
    def _fit(self, x, y, method: Literal['poly', 'interp1d', 'pchip'] = 'pchip', degree: int = 2):
        """
        Fits a curve (e.g., a reservoir elevation-storage curve) using the provided data.

        It supports polynomial fitting, linear interpolation, 
        and shape-preserving cubic Hermite interpolation (PCHIP).
    
        Parameters
        ----------
        x : array_like
            The independent variable data (e.g., elevation values).
        y : array_like
            The dependent variable data (e.g., storage values).
        method : {'poly', 'interp1d', 'pchip'}, optional
            The fitting method to use:
            - 'poly' fits a polynomial of specified degree (e.g., quadratic).
            - 'interp1d' performs linear interpolation.
            - 'pchip' uses shape-preserving cubic Hermite interpolation (default).
        degree : int, optional
            Degree of the polynomial if `method='poly'`. Ignored for other methods. 
            Default is 2.
    
        Returns
        -------
        callable
            A function (or object with a call method) that takes input values (x) 
            and returns estimated output values (y). The return type depends on the method:
            - `np.Polynomial` for polynomial fitting,
            - `scipy.interpolate.interp1d` for linear interpolation,
            - `scipy.interpolate.PchipInterpolator` for PCHIP.
        """
    
        if method.lower() == 'poly':
            coefficients_zv = np.polyfit(x, y, degree)
            curve = np.Polynomial(coefficients_zv[::-1])
        elif method.lower() == 'interp1d':
            curve = interp1d(x=x, y=y, kind='linear', assume_sorted=True)
        elif method.lower() == 'pchip':
            curve = PchipInterpolator(x=x, y=y)
        else:
            raise ValueError(f'"method" must be either "poly", "interp1d" or "pchip": {method} was provided')
        
        return curve
    
    def _lookup_area(self):
        """
        Estimates the area values for the lookup table

        Returns:
        --------
        numpy.ndarray:
            Area values associated to the elevation values in the lookup table
        callabel:
            Preliminary Elevation-Area curve
        """
        # mean values of area and elevation in each bin
        area_avg = (self['storage'].diff() / self['elevation'].diff()).dropna().values
        elev_avg = (self['elevation'].iloc[:-1].values + self['elevation'].iloc[1:].values) / 2

        # fit a curve to the mean area-elevation pairs
        curve_za = self._fit(elev_avg, area_avg, method='pchip')

        # estimate area values for the elevation bins
        area = curve_za(self.elevation)

        return area, curve_za

    def fit(self, method: Literal['poly', 'interp1d', 'pchip'] = 'pchip', degree: int = 2):
        """
        Fits the forward (Elevation -> Storage) and inverse (Storage -> Elevation) 
        curves using the data stored in the ReservoirCurve lookup table.
    
        It supports polynomial fitting, linear interpolation, 
        and shape-preserving cubic Hermite interpolation (PCHIP).
    
        Parameters
        ----------
        method : {'poly', 'interp1d', 'pchip'}, optional
            The fitting method to use:
            - 'poly' fits a polynomial of specified degree (e.g., quadratic).
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
            - `np.Polynomial` for polynomial fitting,
            - `scipy.interpolate.interp1d` for linear interpolation,
            - `scipy.interpolate.PchipInterpolator` for PCHIP.
        
        Raises
        ------
        ValueError
            If an unsupported fitting method is specified.
        """

        if method not in ['poly', 'interp1d', 'pchip']:
            raise ValueError(f'"method" must be either "poly", "interp1d" or "pchip": {method} was provided')
        
        specs = dict(method=method, degree=degree)

        # Elevation-Storage and viceversa
        self.curve_zv = self._fit(self['elevation'], self['storage'], **specs)
        self.curve_vz = self._fit(self['storage'], self['elevation'], **specs)

        if 'area' in self.columns:
            # Elevation-Area and viceversa
            self.curve_za = self._fit(self['elevation'], self['area'], **specs)
            self.curve_az = self._fit(self['area'], self['elevation'], **specs)

            # Area-Storage and viceversa
            self.curve_av = self._fit(self['area'], self['storage'], **specs)
            self.curve_va = self._fit(self['storage'], self['area'], **specs)
        
    def storage_from_elevation(self, elevation: Union[pd.Series, np.ndarray]) -> Union[pd.Series, np.ndarray]:
        """
        Produces a time series of reservoir storage given an elevation time series.
    
        Parameters:
        -----------
        elevation: pandas.Series or numpy.ndarray
            Reservoir elevation data
    
        Returns:
        --------
        storage: pandas.Series or numpy.ndarray
            Estimated reservoir storage data
        """
        return self._convert(
            input=elevation, 
            curve_attr='curve_zv', 
            input_var='elevation', 
            output_var='storage'
            )

    def elevation_from_storage(self, storage: Union[pd.Series, np.ndarray]) -> Union[pd.Series, np.ndarray]:
        """
        Produces a time series of reservoir elevation given a storage time series.
    
        Parameters:
        -----------
        storage: pandas.Series or numpy.ndarray
            Reservoir storage data
    
        Returns:
        --------
        elevation: pandas.Series or numpy.ndarray
            Estimated reservoir elevation data
        """
        return self._convert(
            input=storage,
            curve_attr='curve_vz',
            input_var='storage',
            output_var='elevation'
        )
    
    def area_from_elevation(self, elevation: Union[pd.Series, np.ndarray]) -> Union[pd.Series, np.ndarray]:
        """
        Produces a time series of reservoir area given an elevation time series.
    
        Parameters:
        -----------
        elevation: pandas.Series or numpy.ndarray
            Reservoir elevation data
    
        Returns:
        --------
        area: pandas.Series or numpy.ndarray
            Estimated reservoir area data
        """
        return self._convert(
            input=elevation,
            curve_attr='curve_za',
            input_var='elevation',
            output_var='area'
        )
    
    def elevation_from_area(self, area: Union[pd.Series, np.ndarray]) -> Union[pd.Series, np.ndarray]:
        """
        Produces a time series of reservoir elevation given an area time series.
    
        Parameters:
        -----------
        area: pandas.Series or numpy.ndarray
            Estimated reservoir area data
    
        Returns:
        --------
        elevation: pandas.Series or numpy.ndarray
            Reservoir elevation data
        """
        return self._convert(
            input=area,
            curve_attr='curve_az',
            input_var='area',
            output_var='elevation'
        )
    
    def storage_from_area(self, area: Union[pd.Series, np.ndarray]) -> Union[pd.Series, np.ndarray]:
        """
        Produces a time series of reservoir storage given an area time series.
    
        Parameters:
        -----------
        area: pandas.Series or numpy.ndarray
            Estimated reservoir area data
    
        Returns:
        --------
        storage: pandas.Series or numpy.ndarray
            Reservoir storage data
        """
        return self._convert(
            input=area,
            curve_attr='curve_av',
            input_var='area',
            output_var='storage'
        )
    
    def area_from_storage(self, area: Union[pd.Series, np.ndarray]) -> Union[pd.Series, np.ndarray]:
        """
        Produces a time series of reservoir area given a storage time series.
    
        Parameters:
        -----------
        storage: pandas.Series or numpy.ndarray
            Reservoir storage data
    
        Returns:
        --------
        area: pandas.Series or numpy.ndarray
            Estimated reservoir area data
        """
        return self._convert(
            input=storage,
            curve_attr='curve_va',
            input_var='storage',
            output_var='area'
        )

    def plot(
        self,
        attrs: Optional[pd.Series] = None,
        obs: Optional[pd.DataFrame] = None,
        save: Optional[Union[str, Path]] = None,
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
    
        Parameters:
        -----------
        attrs: pandas.Series (optional)
            An optional Series containing reservoir attributes from GRanD or GDW 
            such as 'DAM_HGT_M', 'ELEV_MASL', 'AREA_SKM' or 'CAP_MCM' to draw 
            reference lines on the plots. Defaults to None.
        obs: pandas.DataFrame (optional)
            An optional DataFrame containing observed reservoir data to overlay on 
            the plots. Defaults to None.
            
        **kwargs: Optional keyword arguments to customize the plot:
            - figsize (tuple, optional): Size of the figure (width, height). 
              Defaults to (10, 10).
    
        Returns:
        --------
        tuple: [plt.Figure, np.ndarray]
            A tuple containing:
            - fig (matplotlib.figure.Figure): The main Matplotlib figure object.
            - axes (np.ndarray): A 2x2 array of Matplotlib axes objects.
        """
        fitted_curves = {
            'z_min': self.z_min,
            'z_max': self.z_max,
            'a_min': self.a_min,
            'a_max': self.a_max,
            'curve_zv': self.curve_zv,
            'curve_za': self.curve_za,
            'curve_av': self.curve_av
        }
        
        fig, axes = plot_reservoir_curves(
            lookup_table=self,
            fitted_curves=fitted_curves,
            attrs=attrs,
            obs=obs,
            **kwargs
        )

        if save:
            plt.savefig(save, dpi=300, bbox_inches='tight')
            plt.close(fig)
        else:
            return fig, axes


def plot_reservoir_curves(
    lookup_table: Optional[pd.DataFrame] = None,
    fitted_curves: Optional[dict] = None,
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
        lookup_table (pd.DataFrame, optional): A DataFrame containing the reservoir 
            characteristic curve data. It must contain the following columns:
            - 'elevation' (float): Reservoir elevation (masl).
            - 'area' (float): Reservoir surface area (km2).
            - 'volume' or 'storage' (float): Reservoir storage volume (hm3).
        fitted_curves: dictionary (optional)
            Curves fitted to the lookup table data
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
    size = kwargs.get('size', 4)
    
    fig, axes = plt.subplots(ncols=2, nrows=2, figsize=figsize, sharex='col', sharey='row')

    if attrs is not None:
        dam_hgt_m, elev_masl, cap_mcm, area_skm = attrs.loc[['DAM_HGT_M', 'ELEV_MASL', 'CAP_MCM', 'AREA_SKM']]

    if fitted_curves:
        z_min, z_max = fitted_curves.get('z_min', 0), fitted_curves.get('z_max', 1)
        a_min, a_max = fitted_curves.get('a_min', 0), fitted_curves.get('a_max', 1)

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
    obs_props = dict(cmap=cmap, s=size, alpha=alpha, zorder=1)
    lookup_props = dict(s=size * 2, lw=.6, c='k', marker='+', alpha=1, zorder=2)
    curve_props = dict(lw=1, c='k', zorder=2)
    
    for j, var_x in enumerate(['elevation', 'area']):
        for i, var_y in enumerate(['storage', 'area']):
            
            ax = axes[i,j]
            if i == 1 and j == 1:
                ax.axis('off')
                continue
                
            # scatter plot of the lookup table
            if lookup_table is not None and all(col in lookup_table.columns for col in [var_x, var_y]):
                label = 'lookup data' if (i == 0 and j == 0) else None
                ax.scatter(lookup_table[var_x], lookup_table[var_y], **lookup_props, label=label)

            # line plot of the fitted curves
            if fitted_curves is not None:
                if var_x == 'elevation':
                    x_values = np.linspace(z_min, z_max, 100)
                    if var_y == 'storage' and fitted_curves.get('curve_zv'):
                        ax.plot(x_values, fitted_curves['curve_zv'](x_values), **curve_props, label='reservoir curve')
                    if var_y == 'area' and fitted_curves.get('curve_za'):
                        ax.plot(x_values, fitted_curves['curve_za'](x_values), **curve_props)
                elif var_x == 'area':
                    x_values = np.linspace(a_min, a_max, 100)
                    if var_y == 'storage' and fitted_curves.get('curve_av'):
                        ax.plot(x_values, fitted_curves['curve_av'](x_values), **curve_props)

            # scatter plot of observed data
            if obs is not None and all(col in obs.columns for col in [var_x, var_y]):
                label = 'observations' if (i == 0 and j == 0) else None
                ax.scatter(obs[var_x], obs[var_y], c=obs.index, **obs_props, label='observations')
                
            for x in var_props[var_x]['ref']:
                ax.axvline(x, **aux_props)
            for y in var_props[var_y]['ref']:
                ax.axhline(y, **aux_props)

            if (i == 1) | (j == 1):
                ax.set_xlabel(var_props[var_x]['label'])
            if j == 0:
                ax.set_ylabel(var_props[var_y]['label'])

    if 'title' in kwargs:
        fig.suptitle(kwargs['title'], y=.93)

    # legend
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='center', bbox_to_anchor=(0.75, 0.25), frameon=False)

    return fig, axes


def save_curve(instance: 'ReservoirCurve', file_path: Union[str, Path]):
    """Saves a fitted ReservoirCurve instance to a pickle file."""
    try:
        if isinstance(file_path, str):
            file_path = Path(file_path)

        # Create directory if it doesn't exist
        file_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(file_path, 'wb') as f:
            pickle.dump(instance, f)
        print(f"Successfully saved ReservoirCurve to {file_path}")
    except Exception as e:
        print(f"Error saving file: {e}")


def load_curve(file_path: Union[str, Path]) -> 'ReservoirCurve':
    """Loads a fitted ReservoirCurve instance from a pickle file."""
    with open(file_path, 'rb') as f:
        instance = pickle.load(f)
    print(f"Successfully loaded ReservoirCurve from {file_path}")
    
    return instance