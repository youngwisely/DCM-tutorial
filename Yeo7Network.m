%% main script for DCM tutorial (for empirical data)

clear; close all; clc

% import Yeo 7 network data
% both hemisphere concatenated
addpath('./data')
BOLD = importdata('Yeo7Net_BOLD_Orig.mat');

% set variables
TR = 0.72;
n = 7;
T = 2400;

% run spectral DCM
%==========================================================================
% set options for the DCM structure
options.nonlinear  = 0;
options.two_state  = 0;
options.stochastic = 0;
options.centre     = 1;
options.induced    = 1;

DCM.options = options;

DCM.a    = ones(n,n);
DCM.b    = zeros(n,n,0);
DCM.c    = zeros(n,1);
DCM.d    = zeros(n,n,0);

DCM.Y.y  = BOLD;
DCM.Y.dt = TR;
DCM.U.u  = zeros(T,1);
DCM.U.dt = TR;

%%
% run spectral DCM
%--------------------------------------------------------------------------
DCM   = spm_dcm_fmri_csd(DCM);
spDCM = DCM.Ep.A; spDCM = spDCM -diag(diag(spDCM));

% make figure
Label = {'Visual';'SomatoMotor';'Dorsal Attention'; 'Ventral Attention'; 'Limbic'; 'Default Mode'; 'Control'};
figure; imagesc(spDCM); title('spectral DCM'); 
set(gca,'Xtick', 1:7, 'XTickLabel',Label, 'XTickLabelRotation',40)
set(gca,'Ytick', 1:7, 'YTickLabel',Label, 'XTickLabelRotation',40)
resolution = 100; cmin = min(spDCM(:)); cmax = max(spDCM(:));
nr_blues    = round(resolution* (abs(cmin) / (abs(cmin) + cmax))); 
nr_reds     = round(resolution* (abs(cmax) / (abs(cmin) + cmax))); 
cMap =[flipud(white2blue(nr_blues)); flipud(red2white(nr_reds))];
colormap(cMap)
xlabel('source')
ylabel('target')

% run regression DCM
Y.y = BOLD;
Y.dt = TR;
rdcm = tapas_rdcm_model_specification(Y,[],[]);
[output,~] = tapas_rdcm_estimate(rdcm,'r',[],1);
rDCM = output.Ep.A; 
rDCM = rDCM - diag(diag(rDCM)); 

% make figure
figure; imagesc(rDCM); title('regression DCM'); 
set(gca,'Xtick', 1:7, 'XTickLabel',Label, 'XTickLabelRotation',40)
set(gca,'Ytick', 1:7, 'YTickLabel',Label, 'XTickLabelRotation',40)
resolution = 100; cmin = min(rDCM(:)); cmax = max(rDCM(:));
nr_blues    = round(resolution* (abs(cmin) / (abs(cmin) + cmax))); 
nr_reds     = round(resolution* (abs(cmax) / (abs(cmin) + cmax))); 
cMap =[flipud(white2blue(nr_blues)); flipud(red2white(nr_reds))];
colormap(cMap)
xlabel('source')
ylabel('target')

%%
% network visualization

nn=1;
for i = 1:7
    source_num = sum(numel(find(rDCM(:,i))));
    if source_num > 0
        if i == 1
            if source_num == 1
                nn = 1;
            else
                nn = 1:source_num;
            end
        else
            nn = nn:nn+source_num-1;
        end

        src(nn) = Label(i);
    
        nn = nn(end)+1;
    end
end

nn=1;
for i = 1:7
    targets = find(rDCM(:,i));
    for j = 1:length(targets)
        tar(nn) = Label(targets(j));
        nn=nn+1;
    end
    
end

w = nonzeros(rDCM(:)');

G = digraph(src,tar,w);
G.Edges.Lwidths = 5*abs(G.Edges.Weight)/max(G.Edges.Weight);


figure;
h = plot(G,'MarkerSize',7,'LineWidth',G.Edges.Lwidths,'Layout','layered');
axis equal;
m = turbo;
colormap(m);

% set(h,'EdgeColor',"#bc0202")
h.EdgeCData = nonzeros(rDCM(:));
set(h,'NodeColor',"#cfcfc4")
set(h,'ArrowSize',10)
% set(he3,'FaceColor',"#2121cf")

box off
set(gca,'Visible','off')
colorbar
