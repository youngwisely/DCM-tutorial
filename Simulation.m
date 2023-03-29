clear;close all; clc
% Simulate timeseries
%==========================================================================
rng('default')

% Define variables and constants
% -------------------------------------------------------------------------
T  = 1200;                             % number of observations (scans)
TR = 1;                               % repetition time or timing
n  = 4;                               % number of regions or nodes
t  = (1:T)*TR;                        % observation times

% set options for the DCM structure
options.nonlinear  = 0;
options.two_state  = 0;
options.stochastic = 0;
options.centre     = 0;
options.induced    = 0;

% set DCM matrices
A   = ones(n,n);
B   = zeros(n,n,0);
C   = zeros(n,n);
D   = zeros(n,n,0);
pP  = spm_dcm_fmri_priors(A,B,C,D,options);


pP.A = randn(n,n)/6;
pP.A = pP.A - diag(diag(pP.A));
pP.C = eye(n,n);
pP.transit = randn(n,1)/16;

% integrate states
% -------------------------------------------------------------------------
U.u  = spm_rand_mar(T,n,1/2)/4;      % endogenous fluctuations
U.dt = TR;
M.f  = 'spm_fx_fmri';
M.x  = sparse(n,5);
x    = spm_int_J(pP,M,U);

% haemodynamic observer
% -------------------------------------------------------------------------
for i = 1:T
    y(i,:) = spm_gx_fmri(spm_unvec(x(i,:),M.x),[],pP)';
end

% observation noise process
% -------------------------------------------------------------------------
e    = spm_rand_mar(T,n,1/2)/8;

% nonlinear system identification (DCM for CSD) over subjects
%==========================================================================
DCM.options = options;

DCM.a    = ones(n,n);
DCM.b    = zeros(n,n,0);
DCM.c    = zeros(n,1);
DCM.d    = zeros(n,n,0);

DCM.Y.y  = y+e;
DCM.Y.dt = TR;
DCM.U.u  = zeros(T,1);
DCM.U.dt = TR;

% provisional inversion
%--------------------------------------------------------------------------
DCM   = spm_dcm_fmri_csd(DCM);
corr(pP.A(:),DCM.Ep.A(:))