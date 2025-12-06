<!-- <img src="./images/logo.png" alt="lisflood-reservoirs logo" width="100" align="center"> -->
<img src="./images/copernicus_logo.png" alt="Logo Copernicus" width="280" align="center"><img src="./images/copernicus_emergency_management.png" alt="Logo CEMS" width="200" align="center">

![Python_3.11](https://img.shields.io/badge/Python-%3E%3D3.11-blue?labelColor=343b41) &nbsp; [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) &nbsp; [![Docs](https://img.shields.io/badge/docs-online-blue.svg)](https://github.com/casadoj/lisflood-reservoirs/wiki)


# LISFLOOD reservoirs

This repository contains tools to run and calibrate different reservoir routines meant to be used in large-scale hydrological models like [LISFLOOD Open Source](https://github.com/ec-jrc/lisflood-code).

Five different reservoir routines are implemented in this repository:

* Linear reservoir (class [`Linear`](./src/lisfloodreservoirs/models/linear.py))
* The routine in the hydrological model [LISFLOOD](https://ec-jrc.github.io/lisflood-model/3_03_optLISFLOOD_reservoirs/) (class [`Lisflood`](./src/lisfloodreservoirs/models/lisflood.py))
* The routine in the hydrological model [CaMa-Flood](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2021MS002944) (class [`Camaflood`](./src/lisfloodreservoirs/models/camaflood.py))
* The routine in the hydrological model [mHM](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2023WR035433) (class  [`mHM`](./src/lisfloodreservoirs/models/mhm.py))
* The reservoir model [Starfit](https://www.sciencedirect.com/science/article/pii/S0022169421008933?via%3Dihub) (class [`Starfit`](./src/lisfloodreservoirs/models/starfit/starfit.py))

Apart from the tools to train and fit these reservoir routines, it contains multiple Jupyter Notebooks to create datasets of reservoir attributes and observed time series in several countries: [US](notebook/ResOpsUS/), [Mexico](./notebook/ResOpsMX/), [Brazil](./notebook/ResOpsBR/), [Spain](./notebook/ResOpsUS/)... These datasets have the same structure as the [CARAVAN](https://github.com/kratzert/Caravan) dataset, and are meant not only as the input data for the reservoir routines in this repository, but also to be used as input for deep learning models.

## Installation

Get a local copy of the repository. You can either download it from GitHub or clone it with Git:

```Bash
git clone https://github.com/casadoj/lisflood-reservoirs.git
```

Move to the root directory of the repository you've just copied:

```Bash
cd <YOUR_PATH>/lisflood-reservoirs/
```

Install the package with PiP:

```Bash
pip install .
```

## Quick start

The repository contains 4 tools to calibrate and run reservoir models. The models included in the repository can be classified in two groups: those that can be calibrated with an iterative process (i.e., a genetic algorithm), and those that are simply fitted using standard `SciPy` tools. To the first group belong the **linear**, **LISFLOOD**, **Camaflood** and **mHM** models; the tools [`run_reservoir`](#run_reservoir) and [`cal_reservoir`](#cal_reservoir) apply to this group. To the second group belongs the **Starfit** model; the tools [`run_starfit`](#run_starfit) and [`fit_starfit`](#fit_starfit) apply to it.

### Configuration

All the tools require a configuration file as the input. A template of this configuration file can be found [here](./src/lisfloodreservoirs/config.yml). The structure of this template is applicable to all the tools.

The configuration file has three sections dedicated to data, simulation, and calibration, repectively.

* The data section defines the location of the reservoir data set and the files that defines the reservoirs to be used (TXT format) and the study period for each of those reservoirs (Pickle format). All the tools are based in a fixed dataset structure:
    *  Attributes must be in a subfolder named *attributes* within the dataset folder.
    *  Time series must be in a subolder named *time_series/csv* within the dataset folder.
* The simulation section defines the reservoir model to be used (`linear`, `lisflood`, `camaflood`, `mhm` or `starfit`) and the folder where the results of the simulation with default parameters will be saved.
* The calibration section defines the name of the input variable, the target or targets of the calibration (`storage`, `outflow` or both), the parameters of the SCE-UA algorithm, and the directory where results will be saved.

### Tools

To run the tools from the command prompt, the instruction is always the same only changing the name of the tools. For instance, to fit the Starfit model:

```Bash
fit_starfit --config-file config.yml
```

#### [`run_reservoir`](./src/lisfloodreservoirs/simulate.py)

This tool simulates the reservoir module with default parameters. It is applicable to the **linear**, **LISFLOOD**, **Camaflood** and **mHM** models.

```
usage: simulate.py [-h] -c CONFIG_FILE [-w]

Run the reservoir routine with default parameters

options:
  -h, --help
                          Show this help message and exit
  -c CONFIG_FILE, --config-file CONFIG_FILE
                          Path to the configuration file
  -w, --overwrite
                          Overwrite existing simulation files. Default: False
```

#### [`cal_reservoir`](./src/lisfloodreservoirs/calibrate.py)

This tool calibrates the reservoir model using the algorithm Shuffle Complex Evolution - University of Arizona (SCE-UA). It is applicable to the **linear**, **LISFLOOD**, **Camaflood** and **mHM** models, and it can calibrate the observed storage, outflow, or both at the same time. Eventually, the model is run with the optimised parameters.

```
usage: calibrate.py [-h] -c CONFIG_FILE [-w]

Run the calibration script with a specified configuration file.
It calibrates the reservoir model parameters of the defined routine using the
SCE-UA (Shuffle Complex Evolution-University of Arizona) algorithm for each of
the selected reservoirs.
The optimal parameters are simulated and plotted, if possible comparing against
a simulation with default parameters

options:
  -h, --help
                          Show this help message and exit
  -c CONFIG_FILE, --config-file CONFIG_FILE
                          Path to the configuration file
  -w, --overwrite
                          Overwrite existing simulation files. Default: False
```

#### [`fit_starfit`](./src/lisfloodreservoirs/fit_starfit.py)

This tool fits the Starfit reservoir model to the observed data.

```
usage: fit_starfit.py [-h] -c CONFIG_FILE [-w]

Fit the storage and release rules for the Starfit reservoir routine.
The fitted models are saved as Pickle files and plotted against the
observed data used for fitting.

options:
  -h, --help
                          Show this help message and exit
  -c CONFIG_FILE, --config-file CONFIG_FILE
                          Path to the configuration file
  -w, --overwrite
                          Overwrite existing model. Default: False
```

#### [`run_starfit`](./src/lisfloodreservoirs/run_starfit.py)

This tool runs the Starfit reservoir model that was previously fitted with the tool [`fit_starfit`](#fit_starfit).

```
usage: run_starfit.py [-h] -c CONFIG_FILE [-w]

Run Starfit simulation with the paremeter fitted using `fit_starfit`.
The simulated time series are saved as CSV files. To analyse the results,
the code creates a CSV file of performance metrics, and a scatter and a 
line plot comparing the observed and simulated time series.

options:
  -h, --help
                          Show this help message and exit
  -c CONFIG_FILE, --config-file CONFIG_FILE
                          Path to the configuration file
  -w, --overwrite
                          Overwrite existing simulation files. Default: False
```

## Datasets

Casado Rodríguez, J., Disperati, J., & Salamon, P. (2025). ResOpsUS+CARS: Reservoir Operations US and CAtchment and Reservoir Static attributes (1.0) [Data set]. European Commission - Joint Research Centre. https://doi.org/10.5281/zenodo.15978041

Casado Rodríguez, J., Disperati, J., & Salamon, P. (2025). ResOpsBR+CARS: Reservoir Operations Brazil and CAtchment and Reservoir Static attributes (1.0) [Data set]. European Commission - Joint Research Centre. https://doi.org/10.5281/zenodo.16096623
