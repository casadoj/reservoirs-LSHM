from setuptools import setup, find_packages

setup(
    name='lisflood-reservoirs',
    version='1.2.0',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    entry_points={
        'console_scripts': [
            'run_reservoir=lisfloodreservoirs.simulate:main',
            'cal_reservoir=lisfloodreservoirs.calibrate:main',
            'fit_starfit=lisfloodreservoirs.fit_starfit:main',
            'run_starfit=lisfloodreservoirs.run_starfit:main',
            'catchstats=lisfloodreservoirs.catchstats:main',
            'ncextract=lisfloodreservoirs.ncextract:main'
        ],
    },
    install_requires=[
        'cartopy',
        # 'cfgrib',
        'dask',
        'gwwapi',
        'lxml',
        'matplotlib',
        'mpi4py',
        'netcdf4',
        'numpy',
        'pandas',
        'pyyaml',
        # 'scikit-learn',
        'seaborn',
        'spotpy',
        'statsmodels',
        'timezonefinder',
        'tqdm',
        'xarray',
    ],
    author='Jesús Casado Rodríguez',
    author_email='jesus.casado-rodriguez@ec.europa.eu',
    description='Package to simulate reservoir operations according to different modelling routines.',
    keywords='hydrology reservoir simulation calibration',
)
