# ResOpsES: Reservoir Operations Spain

ResOps-ES is a hydrometeorological dataset created in this repository that covers 291 reservoirs in Spain and the time period from 1991 to 2023. The dataset includes both reservoir static attributes and time series. The final purpose of ResOpsES is to train hydrological models that reproduce reservoir operations.

The time series were extracted from the [database of the Spanish Ministry of the Environment](https://ceh.cedex.es/anuarioaforos/default.asp), whereas the reservoir characteristics are a combination of the data in the [Inventory of Dams and Reservoirs of Spain](https://www.miteco.gob.es/es/agua/temas/seguridad-de-presas-y-embalses/inventario-presas-y-embalses.html), GRanD and the database of the [International Commission on Large Dams](https://www.icold-cigb.org/).

This document explains the process of creation of the data set _Reservoir Operations Spain_ (ResOpsES) that will be later on used to train the LISFLOOD reservoir routine.

The dataset has been added to this repository in the folder [data/ResOpsES](../../data/ResOpsES/). The structure of the dataset is the following:

* Subfolder _attributes_ contains CSV files with the static attributes. The field _SNCZI_ is the connection between all these files. The number of reservoirs in the CSV files differed depending on the availability of the specific source.

* Subfolder _timeseries_ contains the daily time series associated to the 304 reservoirs selected. The time series are both in CSV or NetCDF format. The file name is the _SNCZI_ code needed to connect the time series and the static attributes. The time series expand from January 1st 1990 (the start of EFASv5 long run) to present, even tough records were available from before. The following time series are included:

   * _level_ (meters above sea level), _storage_ (hm³) and _outflow_ (m³/s) are the oberved time series. Depending on the availability, some of these variables may be missing in some reservoirs.
   * _dis_efas5_ is the inflow time series as specific discharge (mm/d).
   * _e0_emo1_ (mm), _pr_emo1_ (mm) and _ta_emo1_ (°C) are the time series derived from EMO-1: open water evaporation, precipitation and air temperature.
   
> **Note**. Units may need to be changed. For instance, _outflow_ may be converted to mm/d for consistency with _dis_efas_5_.

## 1 Data

The dataset is made up of static attributes that define the characteristics of both the reservoir and its catchment, and time series.

### 1.1 Static attributes

#### 1.1.1 Reservoir attributes

I compared 4 sources of information:

* The [Spanish Inventory of Reservoirs and Dams](https://www.miteco.gob.es/es/cartografia-y-sig/ide/descargas/agua/inventario-presas-embalses.html). This dataset contains mulutiple reservoir and dam attributes for 3170 reservoirs.
* The reservoir attributes in EFASv5 (European Flood Awareness System), which were extracted from GLWD or GRanD. EFASv5 models 245 reservoirs in Spain, but only two attributes of interest (total storage and catchment area).
* [GRanD](https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1890/100125) dataset (Global Reservoirs and Dams). GRanD contains information for 262 reservoirs in Spain, including reservoir storage and area, catchment area, reservoir use, outflow capacity... 
* The ICOLD (Internation Commission on Large Dams) [World Register of Dams](https://www.icold-cigb.org/GB/world_register/world_register_of_dams.asp). ICOLD contains information for 1013 reservoirs in Spain, including reservoir geometry (volume and area), catchment area, reservoir use, etc.

#### 1.1.2 Catchment attributes

To replicate the information in EFASv5, I have created catchment attributes based on:

* [LISFLOOD static maps](https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/CEMS-EFAS/LISFLOOD_static_and_parameter_maps_for_EFAS/). From these maps I computed several statistics (mean, sum, standard deviation, etc) depending on the specific static map. 
* [EMO-1 meteorology](https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/CEMS-EFAS/meteorological_forcings/EMO-1arcmin/). From this dataset I computed monthly and annual averages of precitpitation, air temperature and open water evaporation.

### 1.2 Time series

#### 1.2.1 Reservoir time series

I have found three public sources of **observed** reservoir time series.

* [_Anuario de Aforos de España_](https://www.miteco.gob.es/es/agua/temas/evaluacion-de-los-recursos-hidricos/sistema-informacion-anuario-aforos.html). This is the official national database of reservoir (also gauging stations and canals) observations created by CEDEX (_Centro de Experimentación_). It contains daily time series of level, storage and outflow for **394 reservoirs**.
* [_Agència Catalana de l'Aigua_](https://analisi.transparenciacatalunya.cat/es/Medi-Ambient/Xarxes-de-control-del-medi-consulta-de-l-aigua-i-e/wc95-u57z/about_data). The catchments exclusively inside Catalunya are managed by ACA. The ACA database contains daily time series of percentage filling, level and volume for **15 reservoirs**.
* [_Hidrosur_](http://www.redhidrosurmedioambiente.es/saih/datos/a/la/carta). The catchments exclusively inside Andalucía are managed by the regional government and the data infrastratucture is maintained by Hidrosur. Their website allows access to the historical records, but it's limited to retrievals of 31 days; therefore, I have directly contacted the agency to ask for the reservoir time series. The database contains hourly time series of level and storage for **23 reservoirs**.

All the reservoir time series above were forwarded to the Hydrological Data Colection Centre (HDCC) to be quality checked and converted to UTC+00. At the moment of creation of the ResOpsES dataset, only the CEDEX time series had been added to the HDCC database.

None of the previous sources include reservoir inflow, which is the only dynamical input of the LISFLOOD reservoir routine. For this reason, and to better represent the conditions under which the LISFLOOD reservoir routine works, I have added to the ResOpsES dataset the inflow time series simulated in the EFASv5 long run.

#### 1.2.2 Meteorology

Meterological forcing is, a priori, not needed for the simulation of LISFLOOD reservoir routine. However, I have added daily time series of average precipitation, air temperature and open water evaporation computed from EMO-1.


## 2 Methods

### 2.1 Preprocessing observed time series

The first set of notebooks in this folder were used to preprocess the raw observed time series and prepare them to be sent to HDCC:

* [1.1_CEDEX_reservoirs.ipynb](1.1_CEDEX_reservoirs.ipynb)
* [1.2_Spanish_inventory.ipynb](1.2_Spanish_inventory.ipynb)
* [1.3_ACA_reservoirs.ipynb](1.3_ACA_reservoirs.ipynb)
* [1.4_Hidrosur_reservoirs.ipynb](1.4_Hidrosur_reservoirs.ipynb)

The preprocessing consisted in generating a CSV file with a similar format for each reservoir and a table of attributes for all the reservoirs in a dataset (CEDEX, ACA or Hidrosur) containing some information needed by HDCC to apply the quality checks. 

### 2.2 Selection of reservoirs

The notebook [3.3_select_reservoirs.ipynb](3.3_select_reservoirs.ipynb) combines the four sources of reservoir attributes and the three sources of reservoir time series to define the set of reservoirs that will be included in ResOpsES. This selection of reservoirs was the most time consuming task in the development of the ResOpsES dataset.

As mentioned in the section [Data](#Data), the number of reservoirs in each of the sources varies from more than 3000 in the Spanish Inventory to approximately 250 in EFAS or GRanD. Obviously, many reservoirs are repeated in several databases, but the problem is that there is not a simple way to compare/combine the databases, as they do not share reservoir ID, the coordinates are note exactely the same, and the reservoir/dam name may not be same. In the end, reservoir/dam names proved to be the best approach to connect the different databases; however, making sure that the names in all databases were identical required some manual revision. As CEDEX is the most comprehensive database (over 3000 reservoirs), I have used the CEDEX reservoir ID (called _SNCZI_ in ResOpsES) as the connection betweem databases. In simple terms, I have found the _SNCZI_ code for the reservoirs in EFAS, GRanD and ICOLD.

Once the databases are connected, I could assess the actual number of reservoirs represented across all databases (3272). This number of reservoirs is excessively large to be included in EFAS, and includes many reservoirs that are either very small in terms of total storage, or have a limited effect on the hydrological regime as they control very small catchment areas. Hence, I have filtered the set of reservoirs on two conditions:

* Total storage must be greater than or equal to 10 hm³.
* Catchment area must be greater than or equal to 25 km², which represents approximately 10 pixels in EFASv5.

Both conditions were selected subjectively. The condition on storage proved to be more limiting, especially in the most comprehensive databases such as CEDEX and ICOLD. After filtering, a total of 509 reservoirs were selected. From those, observed time series are available for 304. 

***Table 1**. Summary of the selection of reservoirs for the ResOpsES dataset. Each row represents a different source, and total represents their combination (not the sum). "Initially" indicates the number of reservoirs in the raw database, "filtered" the number of reservoirs that fulfil the storage and catchment area conditions, and "time series" the number of filtered reservoirs for which observed time series are available.*

|           | initially | filtered | time series |
| --------- | --------- | -------- | ----------- |
| CEDEX     |      3170 |      464 |         298 |
| EFAS      |       245 |      242 |         223 |
| ICOLD     |      1013 |      329 |         279 |
| GRanD     |       262 |      253 |         214 |
| **total** |  **3272** |  **509** |     **304** |

In a first version, ResOpsES will include only the 304 reservoirs that fulfil the two conditions and for which records are available, because these are the reservoirs with all the information needed to train the LISFLOOD reservoir routine. In the future, the dataset can be extended to the 509 reservoirs that fulfil the conditions; since records are not available for this reservoirs, they will only be used as PUB (prediction in ungauged basins).

<img src='map_selected_reservoirs.jpg' width='700'>

***Figure 1**. Map of the 304 reservoirs included in the ResOpsES dataset. The size of the dots represent reservoir storage, and the colour the catchment area (yellow for smaller and purple for larger catchments).*

### 2.3 Reservoir attributes

The same notebook [3.3_select_reservoirs.ipynb](3.3_select_reservoirs.ipynb) exports CSV files with the reservoir attributes for each of the data sources. These files are located in the _attributes_ folder of the dataset.

* *attributes_CEDEX.csv*
* *attributes_EFAS.csv*
* *attributes_GRanD.csv*
* *attributes_ICOLD.csv*

Each of these files include the reservoirs in that particular database that fulfil the 2 conditions (column "filtered" in Table 1). The attribute in each file differ, but some common rules apply:

* The field _SNCZI_ is the reservoir code in the Spanish Inventory and is the only connection between all the files.
* Important fields have been renamed for consistency across databases:
    * *RES_NAME* and *DAM_NAME* for reservoir/dam name.
    * *LAT* and *LON*
    * *CATCH_SKM*: catchment area in km².
    * *CAP_MCM*: reservoir storage capacity in hm³ (million cubic meters).
* Reservoir use has been converted into a one hote encoder, i.e., each possible reservoir use is a boolean field. Apart from that, the fields *main_use* and *single_use* identify which of the reservoir uses is the most important (take this with a pinch of salt) or if the reservoir has only one use.

### 2.4 Catchment attributes

The catchment attributes are areal (polygon) statistics of either LISFLOOD static maps or the EMO-1 meteorology.

The first step is to define the catchments of each of the reservoirs included in ResOpsES. After several tries, the most appropriate procedure is simply to apply `cutmaps` on the *upArea.nc* map (upstream area) iteratively on each of the reservoirs. The location of the reservoirs must be coherent with LISFLOOD *ldd.nc* map (local direction drainage); I did this step manually on GIS. The output of `cutmaps` is a set of mask NetCDF files for all the reservoirs in the dataset.

After creating the catchment masks, the notebook [6_create_dataset.ipynb](6_create_dataset.ipynb) generates catchment statistics from the **LISFLOOD static maps**. The statistics are specific for each "category" of static maps: geomorphology, land use, crop coefficient, streams, soil properties, LAI, water demand. They include mean, standard deviation, mininum and maximum, but the functions applied. In the case of time series (LAI and water demand), monthly and yearly aggregates were computed. The result is a table of 110 attributes for the 304 reservoirs saved in the file *attributes_EFAS_static_maps.csv*.

Similarly, I have computed catchment statistics from the **LISFLOOD model parameter maps**. In this case the only statistic is the mean. The result is a table with 14 attributes for the 304 reservoirs saved in the file *attributes_EFAS_model_parameters.csv*.

Moreover, I have added to the ResOpsES dataset catchment statistics from **EMO-1** and **EFAS5**. To generate catchment statistics from EMO-1 I have developed a new tool named `catchstats` in the [lisflood-utilities](https://github.com/ec-jrc/lisflood-utilities) repository. Using this tool I have computed monthly and yearly averages of the three meterological variables in EMO-1: precipitation, air temperature and open water evaporation. From EFAS5, I used the inflow time series extracted from the long run (see following section) to calculate monthly average inflow and average, minimum and maximum annual inflow. In total, there are 54 meterological attributes saved in the file *attribues_EFAS_hydrometeorology.csv*.

### 2.5 Input time series: inflow and meteorology

Reservoir inflow time series were extracted from EFASv5 long run. The tool `ncextract` in the lisflood-utilities repository could be used to carry out this task, but it kept failing at reading some of the EFAS NetCDF files. Instead, as the tool is rather simple, I have applied the same procedure in the notebook [4_EFASv5_inflow.ipynb](4_EFASv5_inflow.ipynb).

The meteorological time series by applying the `catchstats` tool, as first step in the computation of catchment meteorological attributes explained in the previous section.