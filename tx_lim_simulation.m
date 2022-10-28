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

% optional stability test to help speed up the process
option = 0;
if option == 1
    stabilityFlag = real(det((coef)^1e+20));
    if isnan(stabilityFlag) || isinf(stabilityFlag)
        error('The linear system is unstable! Decreasing dt may help!')
    end
end

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

% message for better understanding of the simulation system
if any(real(eig(L))>0)
    disp('Warning: Eigenvalues of L have positive real parts!')
end

if any(real(DD)<0) 
    disp(['SUM over total: ',num2str(sum(real(DD)))])
    disp(['SUM over non-negative: ',num2str(sum(real(DD(ploc))))])
    disp('Warning: Eigenvalues of Q have negative real parts!')
    disp(['Number of non-negative eigenvalues: ', num2str(length(ploc)),'/',num2str(length(DD))])
end
end

function [X0,Xp] = tx_integration(X0,coef,noise,totalMon,dt)

Xp = [];
g = size(X0,2);
for i = 1:totalMon/dt
    tx_disp_progress(i,totalMon/dt,10)
    
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