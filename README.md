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
% Xu et al 2022
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

$\mathbf{LU} = \mathbf{U\Lambda}$ where $\mathbf{U}$ is the matrix of eigenvectors and $\mathbf{\Lambda}$ is the diagonal matrix of eigenvalues, $\lambda_i$. 

$\mathbf{V}$, the eigenvectors of $\mathbf{L}$’s adjoint, is simply determined by $\mathbf{V}^\mathrm{H}=\mathbf{U}^{-1}$, such that $\mathbf{L}^\mathrm{H} \mathbf{V}=\mathbf{VΛ}^\*$, where $\mathrm{H}$ is the conjugate transpose and $*$ is the conjugate. 

Check this [Wolfram World page](https://mathworld.wolfram.com/Eigenvector.html) for useful information to understand the derivation.

The least damped mode corresponds to the mode that decays the slowest, i.e., the magnitude of $|\Re(\lambda_i)|$ is the smallest among all eigenvalues. The associated eigenvector $\mathbf{u}_i$ is the spatial pattern of the trend mode, and the time series of the trend mode is obtained by $\mathbf{v}^{\mathrm{H}}_i \mathbf{x}(t)$.

The trend component in PC space is then determined by multiplying the trend spatial pattern with the temporally varying amplitudes, i.e., $\mathbf{x}_{TR}(t) = \mathbf{u}_i \mathbf{v}^{\mathrm{H}}_i \mathbf{x}(t) $.

After getting $\mathbf{x}_{TR}(t)$ in PC space, we (in Xu et al 2022) obtained the spatial temporal evolution of the trend by dot product of EOFs with $\mathbf{x}\_{TR}(t)$. This spatial temporal evolution is subtracted from original climate variables for detrending. 

The following [MATLAB function](https://github.com/Tongtong-Xu-PSL/LIM/blob/main/tx_lim_trend.m) can be used to extract the LIM trend. 

```Matlab
function [u,alpha,Xtr] = tx_lim_trend(X0,Xtau,tau0)
%-----------------
% Input: X(t), X(t+tau0), and the training lag tau0
% Output: u - spatial pattern of LIM trend mode
%         alpha - temporal evolution of the LIM trend mode
%         Xtr - trend component in PC space.        
%
% Note that the least damped mode should typically be stationary, meaning
% the imaginary part of DD(1) should be 0. Please do not use this code if
% you find DD(1) a complex value.
%
% After getting Xtr in PC space, we (in Xu et al 2022) will get the spatial
% temporal evolution of the trend by dot product of EOFs with Xtr. This
% spatial temporal evolution is subtracted from original climate variables
% for detrending.
%
% Xu et al 2022
%-----------------

% obtain LIM operator
[L,~] = tx_lim_operator(X0,Xtau,tau0);

% eigendecomposition of L matrix: u and adjoint matrix v
[U,D] = eig(L);
V = inv(U)';

% sort eigenvalues according to the decay rate 
% and sort eigenvectors accordingly
DD = diag(D);
[~,loc] = sort(real(DD),'descend');
UU = U(:,loc);
VV = V(:,loc);

% get the least damped mode
u = UU(:,1); % spatial pattern of LIM trend mode
v = VV(:,1);
alpha = v'*X; % time series of LIM trend mode

% trend in PC space
Xtr = u*alpha;

end
```

## 4. Conducting LIM climate simulation

LIM simulations may be generated by integrating white noise forcing with observationally constrained spatial structures (determined from $\mathbf{Q}$) forward in time. Example applications include Capotondi & Sardeshmukh 2017, Xu et al 2021. The following schematic diagram illustrates the process of stochastic integration using LIM,

<img src="https://github.com/Tongtong-Xu-PSL/LIM/blob/main/schematic_simulation_process.png " width="50%" />

Repeating the above process, one can obtain an ensemble of realizations in length of interest. Useful references that document the integration process in detail includes Penland & Matrosova 1994, [supplementary file](https://agupubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1029%2F2020GL090661&file=2020GL090661-sup-0001-Text+SI-S01.pdf) of Xu et al 2021.

The following [MATLAB code](https://github.com/Tongtong-Xu-PSL/LIM/blob/main/tx_lim_simulation.m) can be used to perform ensemble of LIM simulation.

```Matlab
function Xg = tx_lim_simulation(X,tau0,group)
%-----------------
% Input: X(t), the training lag tau0, and the number of ensemble members
% Output: Xg(t), the ensemble with each realization the same length as
%         X(t), and the total number of realizations determined by "group"
%
% Note that "group" needs to be multiple of 20, e.g., 40, 100.
%
% Xu et al 2022
%-----------------

% basic parameters & vectors;
dt = 16/24/30;
samplelen = size(X,2);

% get L & Q from LIM
X0 = X(:,1:end-tau0);
Xtau = X(:,1+tau0:end);
[L,Q] = tx_lim_operator(X0,Xtau,tau0);

% decompose the noise statistics
[V,D] = eig(Q);
[DD,loc] = sort(diag(D),'descend');
V = V(:,loc);

% keep only positive eigenvalues and rescale
ploc = find(real(DD)>=0);
Dp = DD(ploc);
Dp = Dp*sum(real(DD))/sum(real(Dp));
Vp = V(:,ploc);

% matrices used in integration
noise = Vp*diag(sqrt(Dp*dt));
coef = eye(size(L,1))+L*dt;

% initialization (discard the first 2000 years)
subgroup = 20;
X0 = zeros(size(X,1),subgroup);
[X0,~] = tx_integration(X0,coef,noise,12*2000,dt);

% simulation
[~,Xp] = tx_integration(X0,coef,noise,samplelen*group/subgroup,dt);

% organize into an ensemble: each member has length of "samplelen" and in
% total there are "group" number of members
Xg = zeros(size(Xp,1),samplelen,group);
k = 0;
for i = 1:subgroup
    for j = 1:(group/subgroup)
        k = k + 1;
        Xg(:,:,k) = Xp(:,i,(1:samplelen)+samplelen*(j-1));
    end
end

end

function [X0,Xp] = tx_integration(X0,coef,noise,totalMon,dt)

Xp = [];
g = size(X0,2);
for i = 1:totalMon/dt
    % random number generator
    rt = normrnd(0,1,size(noise,2),g);
    
    % integration and iteration
    Xt = coef*X0+noise*rt;
    xp = (X0+Xt)/2;
    
    X0 = Xt;
    if mod(i,1/dt)==0
        Xp = cat(3,Xp,xp);
    end
    if any(isnan(X0))
        disp('Simulation blows up!')
        break
    end
end

end
```

## 5. Final Remarks

The repository provides key elements of LIM analyses in Xu et al 2022: (1) solving LIM operator from observed climate variables, (2) identifying and extracting long-term trend from LIM without assuming whether the trend is linear or nonlinear, (3) conducting large ensembles of climate simulations using LIM. These above codes work under the assumption that the climate system is quasi-stationary and the constructed linear system is stable. For example, [this lecture](http://courses.ece.ubc.ca/491m/lectures/Lecture05.pdf) has some discussion about the stability of linear systems, which could be useful.

Stay tune for more updates!

If you find the materials useful, please help us by citing our very recent paper (Xu et al 2022)!



### References

Penland C, Matrosova L. A Balance Condition for Stochastic Numerical-Models with Application to the El-Nino-Southern Oscillation. Journal of Climate 1994, 7(9): 1352-1372.

Penland C, Sardeshmukh PD. The Optimal-Growth of Tropical Sea-Surface Temperature Anomalies. Journal of Climate 1995, 8(8): 1999-2024.

Alexander M.A, Matrosova L, Penland C, Scott J.D, Chang P. Forecasting Pacific SSTs: Linear inverse model predictions of the PDO. Journal of Climate 2008, 21(2): 385-402.

Frankignoul C, Gastineau G, Kwon YO. Estimation of the SST Response to Anthropogenic and External Forcing and Its Impact on the Atlantic Multidecadal Oscillation and the Pacific Decadal Oscillation. Journal of Climate 2017, 30(24): 9871-9895.

Newman M, Sardeshmukh PD. Are we near the predictability limit of tropical Indo-Pacific sea surface temperatures? Geophysical Research Letters 2017, 44(16): 8520-8529.

Capotondi A, Sardeshmukh PD. Is El Niño really changing? Geophysical Research Letters 2017, 44(16): 8548-8556.

Xu T, Newman M, Capotondi A, Di Lorenzo E. The Continuum of Northeast Pacific Marine Heatwaves and Their Relationship to the Tropical Pacific. Geophysical Research Letters 2021, 48(2): 2020GL090661.

Shin SI, Newman M. Seasonal Predictability of Global and North American Coastal Sea Surface Temperature and Height Anomalies. Geophysical Research Letters 2021, 48(10).

Alexander MA, Shin S-I, Battisti DS. The Influence of the Trend, Basin Interactions, and Ocean Dynamics on Tropical Ocean Prediction. Geophysical Research Letters 2022, 49(3): e2021GL096120.

Xu T, Newman M, Capotondi A, Stevenson S, Di Lorenzo E, Alexander M.A. An increase in marine heatwaves without significant changes in surface ocean temperature variability. Nature Communications 2022. Accepted.

Di Lorenzo E, Xu T, Zhao Y, Newman M, Capotondi A, Stevenson S, et al. Modes and Mechanisms of Pacific Decadal-Scale Variability. Annual Review of Marine Science 2023, 15(1), https://doi.org/10.1146/annurev-marine-040422-084555.