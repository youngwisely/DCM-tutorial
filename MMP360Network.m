%% main script for DCM tutorial (for empirical data)

clear; close all; clc

% import Yeo 7 network data
% both hemisphere concatenated
addpath('./data')
BOLD = importdata('MMP360_BOLD.mat');

% set variables
TR = 0.72;
n = size(BOLD,2);
T = size(BOLD,1);

% run regression DCM
Y.y = BOLD;
Y.dt = TR;
rdcm = tapas_rdcm_model_specification(Y,[],[]);
[output,~] = tapas_rdcm_estimate(rdcm,'r',[],1);
rDCM = output.Ep.A; 
rDCM = rDCM - diag(diag(rDCM)); 

%% or load already analyzed rDCM file
addpath('./results')
rDCM = importdata('rDCM_MMP360.mat');

% make figure
figure; imagesc(rDCM); title('regression DCM'); 
resolution = 100; cmin = min(rDCM(:)); cmax = max(rDCM(:));
nr_blues    = round(resolution* (abs(cmin) / (abs(cmin) + cmax))); 
nr_reds     = round(resolution* (abs(cmax) / (abs(cmin) + cmax))); 
cMap =[flipud(black2cyan(nr_blues)); black2white(nr_reds)];
colormap(cMap)
xlabel('source')
ylabel('target')

outflow = sum(abs(rDCM),1);
inflow = sum(abs(rDCM),2);

surfaceplot(outflow,'MMP','both')
surfaceplot(inflow,'MMP','both')

FC = importdata('FC_MMP360.mat').FC;

outflow = sum(abs(FC),1);
inflow = sum(abs(FC),2);

surfaceplot(outflow,'MMP','both')
surfaceplot(inflow,'MMP','both')

