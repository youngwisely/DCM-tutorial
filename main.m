%% main script for DCM tutorial (for simulation)

clear; close all; clc

%% simulate BOLD signal
% define variables for the simulation
T = 1200;
TR = 0.5;
n = 4;
SNR = 4;

% simulate resting state BOLD signal
[BOLD,GT,options] = Simulation(T,TR,n,SNR);

%% run spectral DCM
% nonlinear system identification (DCM for CSD) over subjects
%==========================================================================
DCM.options = options;

DCM.a    = ones(n,n);
DCM.b    = zeros(n,n,0);
DCM.c    = zeros(n,1);
DCM.d    = zeros(n,n,0);

DCM.Y.y  = BOLD;
DCM.Y.dt = TR;
DCM.U.u  = zeros(T,1);
DCM.U.dt = TR;

DCM   = spm_dcm_fmri_csd(DCM);
spDCM = DCM.Ep.A; spDCM = spDCM -diag(diag(spDCM));


%% run regression DCM
Y.y = BOLD;
Y.dt = TR;
rdcm = tapas_rdcm_model_specification(Y,[],[]);
[output,~] = tapas_rdcm_estimate(rdcm,'r',[],1);
rDCM = output.Ep.A; 
rDCM = rDCM - diag(diag(rDCM)); 

%%
% make a figure
figure; imagesc(GT); colorbar; title('Ground Truth'); 
resolution = 100; cmin = min(GT(:)); cmax = max(GT(:));
nr_blues    = round(resolution* (abs(cmin) / (abs(cmin) + cmax))); 
nr_reds     = round(resolution* (abs(cmax) / (abs(cmin) + cmax))); 
cMap =[flipud(white2blue(nr_blues)); flipud(red2white(nr_reds))];
colormap(cMap)

figure; imagesc(spDCM); colorbar; title('spectral DCM'); 
resolution = 100; cmin = min(spDCM(:)); cmax = max(spDCM(:));
nr_blues    = round(resolution* (abs(cmin) / (abs(cmin) + cmax))); 
nr_reds     = round(resolution* (abs(cmax) / (abs(cmin) + cmax))); 
cMap =[flipud(white2blue(nr_blues)); flipud(red2white(nr_reds))];
colormap(cMap)

figure; imagesc(rDCM); colorbar; title('regression DCM'); 
resolution = 100; cmin = min(rDCM(:)); cmax = max(rDCM(:));
nr_blues    = round(resolution* (abs(cmin) / (abs(cmin) + cmax))); 
nr_reds     = round(resolution* (abs(cmax) / (abs(cmin) + cmax))); 
cMap =[flipud(white2blue(nr_blues)); flipud(red2white(nr_reds))];
colormap(cMap)

% calculate Pearson's correlation coefficient between ground truth and
% posterior expectation of the ground truth
CC1 = corr(GT(:),spDCM(:));
CC2 = corr(GT(:),rDCM(:));

