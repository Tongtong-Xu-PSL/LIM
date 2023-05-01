clear
pp = load('global_sst_ssh_pca_ldm_1958to2017.mat');

sstNum = 29;
sshNum = 22;

%% Set up input for obtaining stochastic simulation result
tau0 = 1;
dt = 16/24/30;
len = length(pp.sst.time);
totalMon = len*3000;

%% stochastic simulation - 1958-2017
X = [pp.sstPCA.PC(1:sstNum,:); pp.sshPCA.PC(1:sshNum,:)];
[Xp,~,~] = tx_lim_stochastic_simulation_noSeason(X,tau0,dt,totalMon,len);

% Xp is in the size of [51, 720, 3000].
% In that 51, 1:29 represents the SST PCs; 30:51 represents the SSH PCs;
% This ensemble can be considered as 3000 members of 60 yrs of global SST & SSH
% simulation (from 60S to 60N).
% Alternatively, if the interest is to have very long simulation, this
% ensemble can also be considered as 20 members of 9000 yrs of global SST &
% SSH simulation (from 60S to 60N); to re-organize this ensemble into the
% size of [51, 9000, 20], it is essentially concatenating Xp(:,:,1) with
% Xp(:,:,2) (i.e., cat(3,Xp(:,:,1),Xp(:,:,2)), and then concatenate the
% result with Xp(:,:,3), so on and so forth, until there is an ensemble of
% 9000 yrs. Then start a new concatenation for the next ensemble.

%% to get the spatial temporal evolution, project the simulation onto EOFs.
% the following uses one ensemble member as the example.
member = 1;
sst2D = pp.sstPCA.EOF(:,1:sstNum)*Xp(1:sstNum,:,member);
ssh2D = pp.sshPCA.EOF(:,1:sshNum)*Xp(sstNum+(1:sshNum),:,member);

% note that the above sst2D & ssh2D is in a matrix of [space, time]; to get
% a matrix in [x,y,time], we need to convert it as followed
sst = ConvertZT_into_XYT(sst2D,pp.sst.mask);
ssh = ConvertZT_into_XYT(ssh2D,pp.ssh.mask);

%% to get the ensemble with trend component, it is simply
sst_with_trend = sst + pp.sst.trend;
ssh_with_trend = ssh + pp.ssh.trend;

% I save Xp so that it is not necessary to run the simulation
save stochastic_simulation_lim_sst29pc_ssh22pc.mat Xp






