%% script to run DCM forward simulation model
clear; close all; clc
% define variables for the simulation
T = 512;
TR = 2;
n = 4;
t  = (1:T)*TR;  

% set options for the DCM structure
options.nonlinear  = 0;
options.two_state  = 0;
options.stochastic = 0;
options.centre     = 1;
options.induced    = 1;

% set DCM matrices
A   = ones(n,n);
B   = zeros(n,n,0);
C   = zeros(n,n);
D   = zeros(n,n,0);
pP  = spm_dcm_fmri_priors(A,B,C,D,options);


A = randn(n,n)/6;
A = A-diag(diag(A));
pP.A = A;
pP.C = eye(n,n);
pP.transit = randn(n,1)/16;
%%
% integrate states
% -------------------------------------------------------------------------
U.u  = spm_rand_mar(T,n,1/2)/4;      % endogenous fluctuations
U.dt = TR;

% plot endogenous fluctuations
figure; plot(U.u); title('endogenous fluctuations');
% plot covariance structure of the endogenous fluctuations
Cov_Endo = corr(U.u);
figure; imagesc(Cov_Endo); title('covariance of endogenous fluctuations')
PSDanalysis(TR,U.u,T);  title('Power spectrum density of endogenous fluctuation'); xlabel('hz'); ylabel('power')
%%
% simulate neuronal state variables
M.f  = 'neuronal_state_function';
M.x  = sparse(n,5);
[x,neural] = Integration(pP,M,U);
PSDanalysis(TR,full(neural)',T); title('Power spectrum density of neuronal signal'); xlabel('hz'); ylabel('power')
%%
y = zeros(T,n);
% haemodynamic observer
% -------------------------------------------------------------------------
for i = 1:T
    y(i,:) = spm_gx_fmri(spm_unvec(x(i,:),M.x),[],pP)';
end
PSDanalysis(TR,y,T); title('Power spectrum density of BOLD signal'); xlabel('hz'); ylabel('power')
% observation noise process
% -------------------------------------------------------------------------
e    = spm_rand_mar(T,n,1/2)/12;

BOLD = y+e;
figure; plot(BOLD); title('simulated BOLD signal')