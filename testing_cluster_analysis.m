% Velocity clustering tests
close all

dataDir   = fullfile('~/Kahuna/data/sabancaya_5_2018/image_exports/');

% 24A
% thermCube   = fullfile(dataDir,'24A/thermCubeAnalysis/thermStats_2020-03-24_z644_x582_t1195.mat');
% velCube     = fullfile(dataDir,'24A/thermCubeAnalysis/opticFlowCNL_20-03-26_n1194.mat');
% idx = [];
% trackParams.trackWindowStart = 1;
% trackParams.trackWindowHeight = 30;
% trackParams.detectionWindowOffset = 18; % (pixels)
% trackParams.detectionWindowHeight = 18; % (pixels)
% trackParams.Tpercentile             = 70;
% trackParams.Gpercetile              = 20;

% 25B 
thermFile   = fullfile(dataDir,'25B/thermCubeAnalysis/thermStats_2020-03-25_z710_x571_t1294.mat');
velCube     = fullfile(dataDir, '25B/thermCubeAnalysis/opticFlowCNL_20-04-09_n1294.mat');
idx = 150;
% trackParams.trackWindowStart = 1;
% trackParams.trackWindowHeight = 50;
% trackParams.detectionWindowOffset = 25; % (pixels)
% trackParams.detectionWindowHeight = 25; % (pixels)
% trackParams.Tpercentile             = 70;
% trackParams.Gpercentile              = 90;
trackParams.trackWindowStart = 110;
trackParams.trackWindowHeight = 55;
trackParams.detectionWindowOffset = 0; % (pixels)
trackParams.detectionWindowHeight = 55; % (pixels)
trackParams.Tpercentile             = 70;
trackParams.Gpercentile              = 90;

use_Tprc = false;
use_Gprc = false;
%% 
% D  = loadif(thermFile,'D');
% V = loadif(velCube,'V');


plotCheckOpticFlow(D,V,idx,trackParams)


Timg = D.T(:,:,idx);
mask = D.mask(:,:,idx);
Vz   = V.Vz(:,:,idx);
Vx   = V.Vx(:,:,idx);
[dTdx,dTdz] = gradient(Timg,D.dx,D.dz);

%% Get roi values
winLims = [trackParams.trackWindowStart trackParams.trackWindowStart+trackParams.trackWindowHeight-1];
trackI = [winLims(1)+trackParams.detectionWindowOffset : ...
                    winLims(1)+trackParams.detectionWindowOffset+trackParams.detectionWindowHeight-1]';
        [roiWin,roiMask,roiPoly] = getROI(mask,'iLims',winLims,'maxRegions',1);
        [roiTrk,trkMask,trkPoly] = getROI(mask,'iLims',trackI,'maxRegions',1);
        
[roi,roimask,roipoly] = getROI(mask,'iLims',trackI,'maxRegions',1);

T = Timg(roimask~=0);
u1 = Vx(roimask~=0);
u2 = Vz(roimask~=0);
G = dTdz(roimask~=0);
%% Percentile filters
Ttop = prctile(T,trackParams.Tpercentile);
Tcut = T(T>=Ttop);
u1c  = u1(T>=Ttop);
u2c  = u2(T>=Ttop);
Gcut = G(T>=Ttop);
Tmask = Timg>=Ttop;

% Could consider trying a distance filter on thermal gradients...

Gpos = G(G>0);
Gtest = [G(G==0); Gpos; -Gpos];
Gthresh = -std(Gtest);

Gmask = dTdz<Gthresh;


%%

%% Plotting
% Plot all initial values in ROI
figure
scatter3(T,G,u2,25)
xlabel('T')
ylabel('G')
zlabel('V')
hold on
scatter3(Tcut,Gcut,u2c,25,[1 0 0])

figure
scatter3(T,u1,u2)
xlabel('T')
ylabel('U')
zlabel('V')
hold on
scatter3(Tcut,u1c,u2c,25,[1 0 0])