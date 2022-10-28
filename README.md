# Linear Inverse Model Simulation

#stochastic-simulation #multivariate-linear-regression

## 1. Introduction

The time evolution of a climate state $\mathbf{x}$ may often be approximated by the stochastically forced linear dynamical system,

$$
\frac{\mathrm{d}\mathbf{x}}{\mathrm{d}t} = \mathbf{Lx} + \mathbf{\xi}
$$

where $\mathbf{x}(t)$ is the climate state, $\mathbf{L}$ is a linear dynamical operator, $\mathbf{\xi}$ is a vector of temporally white noise that may have spatial structure, and $t$ is time. Determining the system from covariances of climate variables results in a Linear Inverse Model (LIM; Penland & Sardeshmukh 1995).

The LIM has been extensively used in predicting seasonal-to-interannual surface ocean conditions (e.g., Alexander et al 2008; Newman & Sardeshmukh 2017; Shin & Newman 2021). But it can also be run as a climate simulation model (Penland & Matrosova 1994). One recent application is Xu et al 2022. And this repository serves as a code source associated with Xu et al 2022.

## 2. Simulation

Some text; test

### References

Penland C, Matrosova L. A Balance Condition for Stochastic Numerical-Models with Application to the El-Nino-Southern Oscillation. Journal of Climate 1994, 7(9): 1352-1372.

Penland C, Sardeshmukh PD. The Optimal-Growth of Tropical Sea-Surface Temperature Anomalies. Journal of Climate 1995, 8(8): 1999-2024.

Alexander M.A, Matrosova L, Penland C, Scott J.D, Chang P. Forecasting Pacific SSTs: Linear inverse model predictions of the PDO. Journal of Climate 2008, 21(2): 385-402.

Newman M, Sardeshmukh PD. Are we near the predictability limit of tropical Indo-Pacific sea surface temperatures? Geophysical Research Letters 2017, 44(16): 8520-8529.

Shin SI, Newman M. Seasonal Predictability of Global and North American Coastal Sea Surface Temperature and Height Anomalies. Geophysical Research Letters 2021, 48(10).

Xu T, Newman M, Capotondi A, Stevenson S, Di Lorenzo E, Alexander M.A. An increase in marine heatwaves without significant changes in surface ocean temperature variability. Nature Communications 2022. In press.