from setuptools import setup, find_packages

setup(
    name='reservoirs-lshm',
    version='1.2.2',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    entry_points={
        'console_scripts': [
            'run_reservoir=reservoirs_lshm.simulate:main',
            'cal_reservoir=reservoirs_lshm.calibrate:main',
            'fit_starfit=reservoirs_lshm.fit_starfit:main',
            'run_starfit=reservoirs_lshm.run_starfit:main',
            'catchstats=reservoirs_lshm.catchstats:main',
            'ncextract=reservoirs_lshm.ncextract:main'
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
        'rioxarray',
        # 'scikit-learn',
        'seaborn',
        'spotpy',
        'statsmodels',
        'timezonefinder',
        'tqdm',
        'xarray',
    ],
    author='Jesús Casado Rodríguez',
    author_email='chus.casado.88@gmail.com',
    description='Package to simulate reservoir operations according to different modelling routines.',
    keywords='hydrology reservoir simulation calibration',
)