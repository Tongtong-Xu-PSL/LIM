function [L,Q] = tx_lim_operator(X0,Xtau,tau0)
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