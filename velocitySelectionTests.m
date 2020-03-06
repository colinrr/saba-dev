%% Velocity selection tests

dataDir   = '~/Kahuna/data/sabancaya_5_2018/image_exports/24A';

dataCube  = fullfile(dataDir, 'thermCubeAnalysis/thermStats_2019-09-18_z641_x591_t1195.mat');

velCube1 = fullfile(dataDir, 'thermCubeAnalysis/velocimetry_20-02-18_n1195_nPyr3_sPyr0-50_nIter3_nSz5_fSz5.mat');
velCube2 = fullfile(dataDir, 'thermCubeAnalysis/opticFlowHS_20-02-20_n1195_nPyr3_sPyr0-50_nIter3_nSz7_fSz15.mat');

idx = 161;

Tprc = 95; % Temperature percentile threshold
V2prc = 70; % Percentile threshold for Horne-Schunk output (velocities are bullshit, but locations look good-ish)
%%
% 
% if ~exist('D','var')
%     disp('Loading data cube...')
%     load(dataCube)
% end
% 
% % if ~exist('V','var')
%     disp('Loading velocity cube 2...')
%     load(velCube2)
%     W2 = -V.Vz(:,:,idx+1);
%     U2 = V.Vx(:,:,idx+1);
% % end
% 
% % if ~exist('V','var')
% disp('Loading velocity cube 1...')
% load(velCube1)
% end

%%

T = D.T(:,:,idx);
imsz = size(Frame);
mask = D.mask(:,:,idx);

poly = mask2poly(mask);

W = -V.Vz(:,:,idx+1);
U = V.Vx(:,:,idx+1);

%% Scaling and thresholding

% Thresholded and normalized Temps
Tsc = (T - vidParams.Tscale(1))./(diff(vidParams.Tscale)).*mask;
Tvals = Tsc(Tsc>0);

Ttop = prctile(Tvals,Tprc);
Vtop = prctile(W,V2prc);

%% 

TscCut = Tsc; TscCut(TscCut<Ttop)=0;
Tmsk = TscCut; Tmsk(Tmsk>0)=1;

Wvals = W(Tmsk==1);

%% 

imagesc(TscCut)
set(gca,'YDir','normal')
