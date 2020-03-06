% test for histogram-tracking

% So far doesn't work bupkis...

clearvars -except D V
close all
%% Testing velocity analysis and selection

dataDir   = '~/Kahuna/data/sabancaya_5_2018/image_exports/24A';

dataCube  = fullfile(dataDir, 'thermCubeAnalysis/thermStats_2019-09-18_z641_x591_t1195.mat');

velCube = fullfile(dataDir, 'thermCubeAnalysis/velocimetry_190918_n1194_nSz5_nPyr3_fSz15.mat');
idx0 = [144]; % Which image to test (can enter multiple to get avg props)

idxTest = [144]; %:5:180];

Tmin = 280;
Tmax = 350; 

% Crude ROI:
roi = [1 75 20 135]; % z1 z2 x1 x2
objectRegion = [55 4 25 25]; % x1 z1 w h

%% Load up and prep data

if ~exist('D','var')
    disp('Loadig data cube...')
    load(dataCube)
end
% if ~exist('V','var')
%     disp('Loadig velocity cube...')
%     load(velCube)
% end

Frame = double(mean(D.T(:,:,idx0),3));
mask = any(D.mask(:,:,idx0),3);
mask = logical(prod(D.mask(:,:,idx0),3));

poly = mask2poly(mask);

% wFrame = sum(V.Vz(:,:,idx+1),3)./numel(idx);

t = D.t;
x = D.x;
z = D.z;

%% Get ROI (window mask, presumably, plus a taper to smooth edges down?)
cutROI = @(im,roi) im(roi(1):roi(2),roi(3):roi(4));

maskROI = cutROI(mask,roi);
% wFrameROI = cutROI(wFrame,roi);
FrameROI = cutROI(Frame,roi);

% wMask = wFrameROI.*maskROI;

xR = x(roi(3):(roi(4)));
zR = z(roi(1):(roi(2)));

% Normalize
maxFrame = max(FrameROI(:));
Tmax = maxFrame;
normFrame = FrameROI;
normFrame(normFrame<Tmin) = Tmin;
normFrame(normFrame>Tmax) = Tmax;
normFrame = (normFrame - Tmin)./(Tmax - Tmin);

% figure
% histogram(wMask(wMask~=0))

%% Plots
% Raw thermal image, ROI
figure('position',[50 200 1300 600])
axa=subplot(1,2,1);
% pcolor(xR,zR,FrameROI)
pcolor(FrameROI)
shading flat
colormap(axa,thermgray(150))
axis equal tight


axa=subplot(1,2,2);
% pcolor(xR,zR,FrameROI)
pcolor(normFrame)
shading flat
colormap(axa,thermgray(150))
axis equal tight
rectangle('Position',objectRegion,'LineWidth',2,'LineStyle','--','EdgeColor',[0.7 0.7 0.7])

%% Do the histogram tracky thing...
tracker = vision.HistogramBasedTracker;
initializeObject(tracker, normFrame , objectRegion);

% I = insertShape(FrameROI,'Rectangle',objectRegion);
% figure
% pcolor(I)

figure
for ii = 1:length(idxTest)
    Framei = double(cutROI(D.T(:,:,idxTest(ii)),roi));
    Tmax = max(Framei(:));
    Framei(Framei<Tmin) = Tmin;
    Framei(Framei>Tmax) = Tmax;
    Framei = (Framei - Tmin)./(Tmax - Tmin);
    
    bbox = step(tracker,Framei);
    
    pcolor(Framei)
    shading flat
    colormap(thermgray(150))
    axis equal tight
    rectangle('Position',bbox,'LineWidth',2,'LineStyle','--','EdgeColor',[0.7 0.7 0.7])
    pause(1)
end
