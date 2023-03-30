function [BOLD,A,options] = Simulation(T,TR,n,SNR)
% Simulate resting state fMRI signal
% Input arguments:
% T = total time points
% TR = repetition time
% n = number of regions
%
% Output
% simulated BOLD signal
%
% Written by Oh Younghyun. Updated March 2023.
%---------------------------------------------------------------------------------------------------------------------
%
% Simulate timeseries
%==========================================================================

% Define variables and constants
% -------------------------------------------------------------------------
t  = (1:T)*TR;                        % observation times
switch SNR
    case 1; sig_x = 4; sig_y = 4; case 2; sig_x = 4; sig_y = 8;
    case 3; sig_x = 4; sig_y = 12; case 4; sig_x = 4; sig_y = 16;
end
        
        

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

% integrate states
% -------------------------------------------------------------------------
U.u  = spm_rand_mar(T,n,1/2)/sig_x;      % endogenous fluctuations
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
e    = spm_rand_mar(T,n,1/2)/sig_y;

BOLD = y+e;



