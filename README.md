# Linear Inverse Model

Some description

## 1. Introduction

The time evolution of a climate state $\mathbf{x}$ may often be approximated by the stochastically forced linear dynamical system,

$$
\frac{\mathrm{d}\mathbf{x}}{\mathrm{d}t} = \mathbf{Lx} + \mathbf{\xi}
$$

where $\mathbf{x}(t)$ is the climate state, $\mathbf{L}$ is a linear dynamical operator, $\mathbf{\xi}$ is a vector of temporally white noise that may have spatial structure, and $t$ is time. Determining the system from covariances of climate variables results in a Linear Inverse Model (LIM; Penland & Sardeshmukh 1995).

The LIM has been extensively used in predicting seasonal-to-interannual surface ocean conditions (e.g., Newman & Sardeshmukh 2017; Shin & Newman 2021). But it can also be run as a climate simulation model (e.g., Xu et al 2022). In this repository, we will detail the process of obtaining climate simulation using LIM, as a code source associated with Xu et al 2022. Please cite Xu et al 2022 if you find the information and codes useful.

## 2. Simulation

Some text

### References

Penland C, Sardeshmukh PD. The Optimal-Growth of Tropical Sea-Surface Temperature Anomalies. Journal of Climate 1995, 8(8): 1999-2024.

Newman M, Sardeshmukh PD. Are we near the predictability limit of tropical Indo-Pacific sea surface temperatures? Geophysical Research Letters 2017, 44(16): 8520-8529.

Shin SI, Newman M. Seasonal Predictability of Global and North American Coastal Sea Surface Temperature and Height Anomalies. Geophysical Research Letters 2021, 48(10).

Xu T, Newman M, Capotondi A, Stevenson S, Di Lorenzo E, Alexander A.M. An increase in marine heatwaves without significant changes in surface ocean temperature variability. Nature Communications 2022. In press.