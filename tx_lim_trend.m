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
% T. Xu
% 2021
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