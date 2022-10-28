# Linear Inverse Model

#stochastic-simulation #multivariate-linear-regression

## 1. Introduction

The time evolution of a climate state $\mathbf{x}$ may often be approximated by the stochastically forced linear dynamical system,

$$
\frac{\mathrm{d}\mathbf{x}}{\mathrm{d}t} = \mathbf{Lx} + \mathbf{\xi}
$$

where $\mathbf{x}(t)$ is the climate state, $\mathbf{L}$ is a linear dynamical operator, $\mathbf{\xi}$ is a vector of temporally white noise that may have spatial structure, and $t$ is time. Determining the system from covariances of climate variables results in a Linear Inverse Model (LIM; Penland & Sardeshmukh 1995).

The LIM has been extensively used in predicting seasonal-to-interannual surface ocean conditions (e.g., Alexander et al 2008; Newman & Sardeshmukh 2017; Shin & Newman 2021). But it can also be run as a climate simulation model (Penland & Matrosova 1994). One recent application is Xu et al 2022. And this repository serves as a code source associated with Xu et al 2022. Codes are written in MATLAB.

## 2. Solving LIM operator

$\mathbf{x}(t)$ usually represents temporally evolving ampitudes of the leading Empirical Orthogonal Functions (EOFs) of climate variables (after removing climatology).

$\mathbf{C}(0) = <\mathbf{x}(t)\mathbf{x}(t)^\mathrm{T}>$ is the auto-covariance matrix of $\mathbf{x}(t)$.

$\mathbf{C}(\tau_0) = <\mathbf{x}(t+\tau_0)\mathbf{x}(t)^\mathrm{T}>$ is the lag-covariance matrix of $\mathbf{x}(t)$. $\tau_0$ is the training lag between $t$ and $t+\tau_0$.

$\mathbf{G}(\tau_0) = \mathbf{C}(\tau_0)/\mathbf{C}(0)$ solves the Green function, which then gives us the LIM operator $\mathbf{L} = \log(\mathbf{G}(\tau_0))$.

The following [MATLAB function](https://github.com/Tongtong-Xu-PSL/LIM/blob/main/tx_lim_operator.m) can be used to solve the LIM operator. Note that it also provides the covariance of the noise $\mathbf{\xi}(t)$, whose role is going to be introduced in the following section.

```Matlab
function [L, Q] = tx_lim_operator(X0,Xtau,tau0)
%-----------------
% Input: X(t), X(t+tau0), and the training lag tau0
% Output: L - LIM operator; Q - noise covariance
% 
% T. Xu
% 2022
%-----------------

C0 = X0*X0'/(size(X0,2)-1);
Ctau = Xtau*X0'/(size(X0,2)-1);

% linear operator
G = Ctau/C0;
L = logm(G)/tau0;

% noise statistics
Q = -1*(L*C0+C0*L');

end
```

## 3. Extracting long-term trend from LIM operator

Several studies (e.g., Frankignoul et al 2017; Alexander et al 2022; Di Lorenzo et al 2023) have shown how the externally forced trend is captured by the least damped eigenmode of $\mathbf{L}$. This is done by performing an eigenanalysis on $\mathbf{L}$; that is,

$\mathbf{LU} = \mathbf{U\Lambda}$ where $\mathbf{U}$ is the matrix of eigenvectors and $\mathbf{\Lambda}$ is the diagonal matrix of eigenvalues ($\lambda_i$). $\mathbf{V}$, the eigenvectors of $\mathbf{L}$’s adjoint, is simply determined by $\mathbf{V}^\mathrm{H}=\mathbf{U}^{-1}$, such that $\mathbf{L}^\mathrm{H} \mathbf{V}=\mathbf{VΛ}^\*$, where $\mathrm{H}$ is the conjugate transpose and $*$ is the conjugate. Check this [Wolfram World page](https://mathworld.wolfram.com/Eigenvector.html) for useful information to understand the derivation.

The least damped mode corresponds to the mode that decays the slowest, i.e., the magnitude of $|\Re(\lambda_i)|$ is the smallest among all eigenvalues. The associated eigenvector $\mathbf{u}_i$ is the spatial pattern of the trend mode, and the time series of the trend mode is obtained by $\mathbf{v}_i^\mathrm{H}\mathbf{x}(t)$. The trend component in PC space is determined by multiplying the sptial pattern with time varying amplitudes, i.e., $\mathbf{x}_{TR}(t)=\mathbf{u}_i\mathbf{v}_i^\mathrm{H}\mathbf{x}(t)$.

### References

Penland C, Matrosova L. A Balance Condition for Stochastic Numerical-Models with Application to the El-Nino-Southern Oscillation. Journal of Climate 1994, 7(9): 1352-1372.

Penland C, Sardeshmukh PD. The Optimal-Growth of Tropical Sea-Surface Temperature Anomalies. Journal of Climate 1995, 8(8): 1999-2024.

Alexander M.A, Matrosova L, Penland C, Scott J.D, Chang P. Forecasting Pacific SSTs: Linear inverse model predictions of the PDO. Journal of Climate 2008, 21(2): 385-402.

Frankignoul C, Gastineau G, Kwon YO. Estimation of the SST Response to Anthropogenic and External Forcing and Its Impact on the Atlantic Multidecadal Oscillation and the Pacific Decadal Oscillation. Journal of Climate 2017, 30(24): 9871-9895.

Newman M, Sardeshmukh PD. Are we near the predictability limit of tropical Indo-Pacific sea surface temperatures? Geophysical Research Letters 2017, 44(16): 8520-8529.

Shin SI, Newman M. Seasonal Predictability of Global and North American Coastal Sea Surface Temperature and Height Anomalies. Geophysical Research Letters 2021, 48(10).

Alexander MA, Shin S-I, Battisti DS. The Influence of the Trend, Basin Interactions, and Ocean Dynamics on Tropical Ocean Prediction. Geophysical Research Letters 2022, 49(3): e2021GL096120.

Xu T, Newman M, Capotondi A, Stevenson S, Di Lorenzo E, Alexander M.A. An increase in marine heatwaves without significant changes in surface ocean temperature variability. Nature Communications 2022. In press.

Lorenzo ED, Xu T, Zhao Y, Newman M, Capotondi A, Stevenson S, et al. Modes and Mechanisms of Pacific Decadal-Scale Variability. Annual Review of Marine Science 2023, 15(1): null.