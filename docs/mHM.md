# The reservoir routine in mHM
*** 

## Reservoir release

The reservoir release at day $t$ is partitioned in a component based on the hedged demand ($\hat{D}_t$) and a component based on the inflow ($I_t$):

$$Q_t = \rho \cdot \kappa_t \cdot \hat{D}_t + (1 - \rho) \cdot I_t$$

The partition coefficient ($\rho$) depends on the degree of regulation of the reservoir ($DOR$), i.e., the ratio between the reservoir storage capacity and the annual inflow. The release of highly regulated reservoirs ($DOR \geq DOR*$) depends only on the demand component ($\rho = 1$), whereas in the rest of the cases the release depends on both the demand and the current inflow components:

$$\rho = min\left(1, \left( \frac{DOR}{DOR*} \right)^\beta \right)$$

$DOR*$ and $\beta$ are model parameters in mHM. Shin et al. (2019) and Hanasaki et al. (2006) used fixed values as reported in Table 1.

In the release formulation, the current filling ($V_t$) of the reservoir is used indirectly to limit the supplied demand by the time varying coefficient $\kappa_t$:

$$\kappa_t = \left( \frac{V_t}{\gamma \cdot V} \right) ^\lambda = \left( \frac{V_t}{V_n} \right) ^\lambda$$

where $V_n$ is the normal filling of the reservoir (mHM estimates this value with the parameter $\gamma$), and $\lambda$ is a parameter that further controls demand hedging. Again, Shin et al. (2019) and Hanasaki et al. (2006) used fixed values for $\gamma$ and $\lambda$ (see Table 1).

## Demand hedging

The hedged demand ($\hat{D}_t$) is a transformation of the actual daily demand ($D_t$) based on the average inflow ($\overline{I}$) and average demand ($\overline{D}$). The $\omega$ parameter represents the annual excess of water with repect to demands and it is used to classify reservoirs in water-stressed or non-water-stressed.

In reservoirs with high water stress ($\frac{\overline{D}}{\overline{I}} \gt 1 - \omega$), the hedged demand is a combination of a fixed percentage of the mean inflow, and a varying term depending on the ratio between the current and the average demand.

$$\hat{D}_t = \omega \cdot \overline{I} + \left( 1 - \omega \right) \frac{D_t}{\overline{D}} \cdot \overline{I}$$

In reservoirs with lower water stress ($\frac{\overline{D}}{\overline{I}} \lt 1 - \omega$), the hedged demand is the average surplus of inflow plus the current demand. 

$$\hat{D}_t = \overline{I} - \overline{D} + D_t$$

mHM considers $\omega$ a model parameter, but previous studies used fixed values (see Table 1).

> **Note**. The demand is further hedged in the computation of the reservoir release depending on the current filling of the reservoirs ($\kappa_t$)

## Model parameters

***Table 1**. Model parameters in the reservoir routine in mHM and values used in previous studies.*

| parameter | description | range | Hanasaki (2006) | Shin (2019) | Sadki (2023) | Shrestha (2024) |
| --------- | ----------- | ----- | --------------- | ----------- | ------------ | --------------- |
| $\omega$  | Divides reservoir according to water stress and is used in demand hedging of the stressed reservoirs | [0, 1] | 0.5 | 0.1 | calibrated | calibrated |
| $DOR*$ | Threshold on the degree of regulation | (0, inf) | 0.5 | 1.176 | calibrated | calibrated |
| $\beta$ | Indirectly controls the partition of releases in demand-based and inflow-based | | 2 | 1 | calibrated | calibrated |
| $\gamma$ | Fraction filled representing the normal storage | [0, 1] | - | 0.85 | calibrated | calibrated |
| $\lambda$ | Controls hedging based on current reservoir filling | | 1 | 1 | 1 | calibrated |

> **Note**. The paper mentions that the mHM reservoir routine has 8 parameters, but from the formulation I only identify 5. It could be that they consider parameters the storage capacity ($V$), average inflow ($\overline{I}$) and average demand ($\overline{D}$), which need to be defined when instantiating a reservoir, but those are not calibration parameters.