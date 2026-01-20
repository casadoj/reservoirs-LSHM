# [Starfit](https://github.com/IMMM-SFA/starfit)


## Introduction

I have implemented in Python the reservoir model introduced in [Turner et al. (2021)](https://www.sciencedirect.com/science/article/pii/S0022169421008933). The paper refers to an R repository called [Starfit](https://github.com/IMMM-SFA/starfit) that includes the functions to fit the storage and release functions of this reservoir routine. I have translated these functions into Python and created a class [`Starfit`](https://github.com/casadoj/lisflood-reservoirs/blob/ResOpsUS/src/lisfloodreservoirs/models/starfit/starfit.py) that takes the fitted parameters and simulates the reservoir.

The full implementation of the Starfit reservoir routine requires the following variables:

* Reservoir attributes:
    * $S$: total reservoir storage [hm3]
    * $\bar{I}$: average reservoir inflow [hm3/week]
* Daily time series:
    * $S_t$: reservoir storage [hm3]
    * $I_t$: reservoir inflow [m3/s]
    * $R_t$: reservoir release [m3/s]
    
The model is based on fitting harmonic functions to the observed storage and outflow time series in order to learn the seasonal behaviour of the reservoir. To avoid gaps and errors in daily data, the fitting is done on weekly aggregated time series (all in hm3 units) and gaps are tried to be filled in by mass balance. The fact that the model is fitted on weekly data does not prevent its use on daily simulations. As we will see later, the harmonic functions include a frequency term ($\omega$); by changing this term we can use the same fitted functions for any temporal resolution.

The functions are fitted on standardised variables to allow regionalization of the parameters. Storage is converted into fraction filled by dividing it by the total storage capacity ($S$), and inflow/release are normalised by the mean inflow ($\bar{I}$).

$$
\begin{align}
S_{st,t} &= \frac{S_t}{S} \\
R_{st,t} &= \frac{R_t - \bar{I}}{\bar{I}} \\
I_{st,t} &= \frac{I_t - \bar{I}}{\bar{I}}
\end{align}
$$

## Storage normal operating range (NOR)

The model uses two harmonic functions to define the normal operating range (NOR) of the standardised storage ($\hat{S}$). These two rules define the upper and lower bounds of the normal reservoir filling for every week of the year. During a simulation, the routine will try to keep the reservoir in the NOR by using a different release operation depending on the storage zone (below, within or above NOR).

$$
\begin{align}
NOR_{up} &= \min \left( \max \left( A + B \cdot \sin 2 \pi \omega t + C \cdot \cos 2 \pi \omega t, S_{st,min} \right) , S_{st,max} \right) \\
NOR_{low} &= \min \left( \max \left( a + b \cdot \sin 2 \pi \omega t + c \cdot \cos 2 \pi \omega t, s_{st,min} \right) , s_{st,max} \right)
\end{align}
$$

Each of the NOR harmonics has 5 parameters: 3 defining the harmonic ($A$, $B$, $C$ in the upper bound), and 2 capping the maximum and minimum values of that NOR ($S_{st,max}$, $S_{st,min}$ in the upper bound). Therefore, the storage routine has 10 parameters.

> $\omega$ is the frequency. Fitting (weekly): $\omega=\frac{1}{52}$. Simulation (daily): $\omega=\frac{1}{365}$<br>
> $t$ is the time of the year. Fitting (weekly): epistemologic week of the year. Simulation (daily): day of the year

To fit these functions only the 3 extreme values for each week of the year are used. Figure 1 clarifies the procedure. It shows the observed reservoir filling for each of the 52 weeks of the year. The highest 3 values for each week (green dots) are used to fit the upper bound of the NOR, and the lowest 3 values (red dots) to fit the lower bound of the NOR.

<img src="../results/ResOpsUS/starfit/355_NOR.jpg" width="600">

***Figure 1**. Fitted normal operating range (NOR) for reservoir 355.*

## Release function

The release function operates differently depending on the storage zone:


$$
R_t = \begin{cases}
R_{min} & if \quad S_{st,t} < NOR_{low} \\
\min \left( \bar{I} \cdot R_{st,t} + \bar{I} , R_{max} \right) & if \quad NOR_{low} \leq S_{st,t} \leq NOR_{up} \\
\min \left( S \cdot \left( S_{st,t} - NOR_{up} \right) + I_t , R_{max} \right) & if \quad S_{st,t} > NOR_{up}
\end{cases}
$$

* **Below NOR**, the release is constant ($R_{min}$). In the original implementation, this minimum release is the 5% quantile of daily inflows. I have changed this definition to the 1% quantile because the original definition did not produce good results in simulation mode.

* When the reservoir filling is **within the NOR**, the routine models the standardised release ($R_{st,t}$) as the sum of an harmonic model ($R_{harm,t}$) that defines the seasonality, and a linear model ($\epsilon_t$) that applies a correction based on the current reservoir filling ($S_{st,t}$) and normalised inflow ($I_{st,t}$).

$$
\begin{align}
R_{st,t} &= R_{harm,t} + \epsilon_t \\
R_{harm,t} &= d \cdot \sin 2 \pi \omega t + e \cdot \cos 2 \pi \omega t + f \cdot \sin 4 \pi \omega t + g \cdot \cos 4 \pi \omega t \\
\epsilon_t &= h + i \cdot \frac{S_{st,t} - NOR_{low}}{NOR_{up} - NOR_{low}} + j \cdot I_{st,t}
\end{align}
$$

> The release harmonic function allows for two periods of higher release throughout the year, whereas the storage harmonic only allows for one period of higher reservoir filling.

> The linear model is fitted on the residuals of the harmonic model. If the $R^2$ of the fitted model is below a threshold (0.2 in the code, 0.3 in the paper), the linear model is discarded; this happens in many reservoirs.

> The actual release is capped by a maximum ($R_{max}$). In the original implementation, this maximum release is the 95% quantile of daily inflows. I have changed this definition to the 99% quantile because the previous value did not work well in simulation mode.

* **Above NOR**, the release tries to bring back the reservor to a normal status by releasing the sum of the inflow and the excess volume.

The above calculated release ($R_t$ ) needs to be checked in terms of mass balance:

$$R_t = \max \left( \min \left( R_t, I_t + S_t \right), I_t + S_t - S \right)$$

In total, the release model has 9 parameters: 4 related to the harmonic function ($d$, $e$, $f$, $g$), 3 to the linear model ($h$, $i$, $j$), and the minimum and maximum releases ($R_{min}$, $R_{max}$). Unlike the storage model, the release model is fitted using the complete set of available weekly releases. Figure 2 shows the harmonic and linear functions of release fitted for reservoir 355.

<img src="../results/ResOpsUS/starfit/355_release.jpg" width="900">

***Figure 2**. Fitted release model for reservoir 355. On the left panel, the harmonic, minimum and maximum releases are shown against weekly observations. On the right panel, the linear model of release residuals based on standardised inflow (X axis) and storage availability (colour lines).*

### Modifications

#### Release within NOR

In several reservoirs, the minimum release ($R_{min}$) is above the minimum of the harmonic release ($R_{harm}$). This means that, at some points of the year, the release within NOR is lower than the release below NOR, which is nonesense.

I have added a new constraint to the release within NOR, to be at least equal to the minimum release:

$$R_{NOR} = \max \left( \min \left( \bar{I} \cdot \left( R_{harm,t} + \epsilon_t \right) + \bar{I}, \, R_{max} \right) , \, R_{min} \right)$$

#### Release below NOR

The original definition generates very noisy releases when storage fluctuates around the lower bound of the normal operating rule ($NOR_{low}$). If storage is slightly over $NOR_{low}$, the release is a function of the harmonic and the linear models; however, when the storage is slightly below $NOR_{low}$, the release is a constant $R_{min}$. It means that from one step to the next there is a big difference in release.
    
I have changed the release function below $NOR_{low}$ to be the minimum value between the inflow ($I_t$) and a linear function of the reservoir filling ($S_{st,t}$). In this way, the changes in release are smooth, and the release is at most equal to the inflow, so the reservoir does not keep emptying:
    
$$R_{low} = \min \left( R_{min} + \left( R_{NOR} - R_{min} \right) \cdot \frac{S_{st,t}}{NOR_{low}} , \, I_t \right)$$

#### Release above NOR

The original routine tries to bring back the reservoir to the upper bound of NOR in one time step, if that release does not exceed the maximum ($R_{max}$). That behaviour often causes noisy releases when the value of the maximum release is relatively large. Following the reasoning used before, I have changed the definition of the release in the flood zone to a linear function of the reservoir filling ($S_{st,t}$).

$$R_{up} = R_{NOR} + \left( R_{max} - R_{NOR} \right) \cdot \frac{S_{st,t} - NOR_{up}}{1 - NOR_{up}}$$

## Results

### Model parameters

#### Storage model

The results of fitting the storage model are available at _Z:/nahaUsers/casadje/datasets/ResOpsUS/results/starfit/NOR/_. For each reservoir, there's a Pickle file with the fitted parameters, and an image similar to that shown in Figure 1.

#### Release model

The results of fitting the release model are available at _Z:/nahaUsers/casadje/datasets/ResOpsUS/results/starfit/release/_. For each reservoir, there's a Pickle file with the fitted parameters, and an image similar to that shown in Figure 2.

### Simulation with fitted parameters

I have run the Starfit model using the fitted parameters and two input data: the observed inflow daily time series, and the observed storage at the beginning of the simulation. 

<img src="../results/ResOpsUS/starfit/355_line.jpg" width="900">

***Figure 3**. Simulation of reservoir 355 with the fitted Starfit model compared with observations.*

This run is comparable to the bivariate calibration of the other reservoir models, as fitting Starfit requires both storage and outflow records. There is a fundamental difference, though, between calibrating the other models and fitting Starfit. The calibration of the other models was done with a genetic algorithm, which implies multiple runs and a relatively long calibration time (several hours in the HPC). Fitting Starfit uses ordinary least squares, with takes seconds in a regular PC. 

Figure 4 shows the comparison of the performance of Starfit with respect to the other 4 reservoirs models included in `lisflood-reservoirs`.

<img src="../results/ResOpsUS/ecdf_KGE_SCEUA-QS.jpg" width="900">

***Figure 3**. Simulation of reservoir 355 with the fitted Starfit model, compared with observations.*

Starfit ranks second, with a performance very similar to Hanazaki. Hanazaki performns slightly better in terms of storage, but slightly worse in terms of outflow. mHM is still the best performing model.

## Extra

The `Starfit` repository includes a function to regionalize model parameters to reservoirs lacking data based on proximity and reservoir use. I haven't converted this function to Python yet.