# [ResOpsUS](https://www.nature.com/articles/s41597-022-01134-7): Reservoir Operations USA

This folder contains several analysis done using as observations the ResOpsUS data, a dataset that includes records for 679 major reservoirs across the US. The time series include inflow, storage, outflow and evaporation, although not all variables are available for all reservoirs. 

The folder [EDA](./EDA) contains the exploratory data analysis of the data in the ResOpsUS dataset and the reservoir simulations in GloFASv4.

The folder [dataset](./dataset) contains the notebooks that creates a new version of the ResOpsUS dataset that includes, apart from the original observations, static attributes from GRanD, GloFAS and ERA5, and time series from GloFAS and ERA5. This new version of the dataset will be used to train and compare reservoir models.

The folder [models](./models) contains notebooks used for the development of the different reservoir routines contained in the repository, and a notebook that compares the performance of these routines.