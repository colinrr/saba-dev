

%% TEST INPUT

% ##############################################################
% Using calculated velocities from Optic Flow to track features
% ##############################################################

% clear all; close all
clearvars -except D V vidParams; %close all
clear textprogressbar
%%
% cubeDir  = '~/Kahuna/data/sabancaya_5_2018/image_exports/24A/thermCubeAnalysis/';
% cubeFile = fullfile(cubeDir,'thermStats_2019-09-18_z641_x591_t1195.mat');
% velCube  = fullfile(cubeDir,'velocimetry_190918_n1194_nSz5_nPyr3_fSz15.mat');
% velCube = fullfile(cubeDir, 'velocimetry_20-02-20_n1195_nPyr3_sPyr0-50_nIter3_nSz5_fSz10.mat');
% velCube = fullfile(cubeDir, 'velocimetry_20-02-18_n1195_nPyr3_sPyr0-50_nIter3_nSz3_fSz15.mat');
% velCube = fullfile(cubeDir, 'velocimetry_20-02-18_n1195_nPyr3_sPyr0-50_nIter3_nSz3_fSz15.mat');
% velCube = fullfile(cubeDir, 'opticFlowFB_20-02-26_n1195_nPyr3_sPyr0-50_nIter3_nSz3_fSz7.mat');
% velCube = fullfile(cubeDir, 'opticFlowFB_20-03-01_n1195_nPyr3_sPyr0-50_nIter3_nSz3_fSz10.mat');
% velCube = fullfile(cubeDir, 'opticFlowFB_20-03-01_n1195_nPyr3_sPyr0-50_nIter3_nSz5_fSz15.mat');

muWinHeight = 30; % Height of windows for mean image profile

winHeight  = 30; % Height of the MAIN window to retrieve thermal stats (px)
winDur     = 1; % Number of frames to average over
tOver      = 0 ; % Number of frames overlapping in time window
srcIdxM    = 1; % Height index of first window (px)

% ---- Source pulse detection parameters (thermPulseDetection.m) ----

% Because velocities and detections are at the leading edge of a moving
% pulse, but thermal statistics should be retrieved from the pulse body,
% the detection window should be vertically offset above the main window.
% 2 parameters below are in units of pixels. Defaults are:
% (ie detection_window_offset  = round(winHeight * 0.5) - tracking starts half a window height up) 
% (ie detection_window_scale  = round(winHeight * 0.5)  - tracking window is one half the main window) 
% detection_window_offset = 23; % Vertical shift of detection window (pix)
% detection_window_scaple  = 18; % Vertical scale of detection window (pix)

det_field = 'prctile'; % Which field of T0 (source window thermal stats) to use for detection
det_chan  = 5; % ONLY if det_field == 'prctile'. Selects which data row to use

% Source signal pre-filter parameters
prefilt_params.taperLength = [];

% STA/LTA detection parameters
det_params.l_sta         = []; %round(det_params(1)*Fs);
det_params.l_lta         = []; %round(det_params(2)*Fs);
det_params.th_on         = []; %det_params(3);
det_params.th_off        = []; %det_params(4);
det_params.min_dur       = []; %det_params(5);
det_params.lta_mode      = [];

trackParams.Tpercentile             = 70;
trackParams.Gpercetile              = 20;
trackParams.detection_window_offset = 18; % (pixels)
trackParams.detection_window_height = 18; % (pixels)
% -------------------------------------------------------------------
% figdir = '/Users/crrowell/Kahuna/data/sabancaya_5_2018/image_exports/24A/figures';


trkIdxM    = 2;
% srcROI = []; % Maybe something similar...

vI0 =2; % Which velocity track(s) to follow?


detection_plotflag  = true; % Plot detection results
filter_plotflag     = true; % Plot pre- and post-filter velocity w/ spectra
% Plot colors
fs = 14;
fc = [0 0 0];

Tedges = linspace(200,420,50);

atmofile = '~/Kahuna/data/sabancaya_5_2018/MODIS_atmospheric_profiles/atmo.mat';
%% DO THE THING
fprintf('\n========= Thermal Pulse Tracking =========\n')
% disp('Loading dat...')

% Load if not already in workspace

if ~exist('V','var')
    disp('Loading velocity data cube...')
    load(velCube)
end
if ~exist('D','var')
    disp('Loading thermal data cube...')
    try 
        load(V.dataCube)
    catch
        error('Thermal data cube associated with this velocity file is missing!')
    end
end
load(atmofile)

%% !!!!! Quick and dirty fix - should actually apply this to the velocimetry !!!
% if size(V.Vz,3)~=numel(D.t)
%     V.Vx = cat(3,zeros([size(V.Vx(:,:,1)) numel(D.t)-size(V.Vz,3)]),V.Vx);
%     V.Vz = cat(3,zeros([size(V.Vz(:,:,1)) numel(D.t)-size(V.Vz,3)]),V.Vz);
% end
V.t = D.t(1:end-1) + diff(D.t)/2;

% !!!!!!
%%
M = size(D.T,1);
N = size(D.T,2);
P = size(D.T,3);

dt = mean(diff(D.t));

Tbins = Tedges(1:end-1) + diff(Tedges)/2;

% Get atmospheric profile relative to vent height
Mz0 = 5243; % Camera altitude
ventz = 5919;

atmoZ = atmo.Height(~isnan(atmo.Height))-ventz;
%% Source function - thermal max, variance, mean/rms/median?
disp('Getting source and detection window statistics...')

% Get time windows
[src_tI,winTime] = getSTFTColumns(P,winDur,tOver,1/dt); % Something weird with this wintTime
if winDur>1
    winTime = mean(D.t(src_tI),1);
else
    winTime  = D.t';
end
Pwin             = size(src_tI,2); % Number of time windows

src_zI = repmat((srcIdxM:srcIdxM+winHeight-1)', [1 Pwin]); % Set source region
det_zI1 = srcIdxM+trackParams.detection_window_offset-1; % Bottom height of detection window
det_zI = repmat((det_zI1+1:det_zI1+trackParams.detection_window_height)', [1 Pwin]); % Set detection window

% trk_zI = repmat((trkIdxM:trkIdxM+winHeight-1)', [1 Pwin]); % Set tracking source region

% Get source window stats
T0 = getThermStats(D.T,D.mask,D.z,D.t,src_zI,src_tI);
V0 = getThermStats(V.Vz,D.mask(:,:,2:end),D.z,V.t,src_zI(:,2:end),src_tI(1:end-1));
% [V0f,V0spec] = filterVelocities(V0,dt,'low',3,4,true);

% Detection window stats
Td = getThermStats(D.T,D.mask,D.z,D.t,det_zI,src_tI);
Vd = getThermStats(V.Vz,D.mask(:,:,2:end),D.z,V.t,det_zI(:,2:end),src_tI(1:end-1));
[Vdf,Vdspec] = filterVelocities(Vd,dt,'low',2,4,filter_plotflag);
% Tmax = squeeze(max(D.T.*D.mask,[],2));


%% Implement detection here
disp('Running source pulse detection...')
if strcmp(det_field,'prctile')
    Tdetect = Td.(det_field)(det_chan,:);
else
    Tdetect = Td.(det_field);
end
[tTrig,yTrig] = thermPulseDetection(Tdetect, Td.t, prefilt_params, det_params, detection_plotflag);

% Get detection start frames
[~,tI0] = closest(tTrig(:,1),D.t);
facs = factor(ceil(numel(tI0)/4)*4);
dims = [prod(facs(1:ceil(numel(facs)/2))) prod(facs(ceil(numel(facs)/2)+1:end))];
rows = min(dims); cols = max(dims);

% Plot Check Detection Frames
% plotROI = zeros(size(tTrig,1),4);
% for pr = 1:size(plotROI,1)
%     plotROI(pr,1) = max([find( sum(D.mask(src_zI(:,1),:,tI0(pr)),2),1,'first' )  src_zI(1,1)  ]);
%     plotROI(pr,2) = min([find( sum(D.mask(src_zI(:,1),:,tI0(pr)),2),1,'last' )   src_zI(end,1)  ]);
%     plotROI(pr,3) = find( sum(D.mask(src_zI(:,1),:,tI0(pr)),1),1,'first' );
%     plotROI(pr,4) = find( sum(D.mask(src_zI(:,1),:,tI0(pr)),1),1,'last' );
% end
% roiX1 = find( sum(sum(D.mask(src_zI(:,1),:,:),3),1),1,'first' );
% roiX2 = find( sum(sum(D.mask(src_zI(:,1),:,:),3),1),1,'last' );
% plotROI = [src_zI(1,1) src_zI(end,1) roiX1 roiX2];

roiX1 = find( sum(sum(D.mask(det_zI(:,1),:,:),3),1),1,'first' );
roiX2 = find( sum(sum(D.mask(det_zI(:,1),:,:),3),1),1,'last' );
plotROI = [det_zI(1,1) det_zI(end,1) roiX1 roiX2];

% plotCheckOpticFlow(D,V,tI0(2:6),plotROI)


%% Use velocities to track thermal features
disp('Retrieving velocity tracks...')
% winI = [33:5:50 123:5:138 208:10:248 286:10:336 394];
% winI = [33:5:58  78:20:500];
% winI = [48 145 215 310 415 498 ];
% winI = [48 145 215 310 415 498 ];

% Vpos = trackVelocities(V,D,src_zI(:,1),tI0(1),trackParams); 
Vpos = trackVelocities(V,D,src_zI(:,1),tI0(1:2),trackParams); 

%% Get stats of tracked thermal features
disp('Retrieving thermal profiles...')

for kk=length(Vpos) %numel(winI):-1:1
    T(kk) = getThermStats(D.T,D.mask,D.z,D.t,Vpos(kk).zI,Vpos(kk).tI);
end

%% Get time-averaged image
muI = [123 500];
Dmu = D;
Dmu.t = Dmu.t(muI(1):muI(2));
Dmu.T = mean(Dmu.T(:,:,muI(1):muI(2)),3);
Dmu.mask = mean(Dmu.mask(:,:,muI(1):muI(2)),3)>0.5;

% Crude but meh
% zImu = zeros(size(src_zI,1),250); % Manual 250, about pixel height of mean plume mask
% zImu(:,1) = src_zI(:,1);
zImu = repmat([1:muWinHeight]',1,250);
zadd = repmat(0:249,size(zImu,1),1);
zImu = zImu + zadd;

Tmu = getThermStats(Dmu.T,Dmu.mask,Dmu.z,Dmu.t(1),zImu,ones(1,size(zImu,2)));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------ (nefarious) PLOTTING DOWN HERE -----------------------
% vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
%% Source history plot

% tStart = [D.t(winI)'; D.t(winI)'];
% % Trange = repmat([min(T0.prctile(:)); max(T0.prctile(:))], [1 size(tStart,2)] );
% Trange = repmat([min(T0.var(:)); max(T0.var(:))], [1 size(tStart,2)] );

% SOURCE WINDOW
sax=plotThermStats(T0,V0,tTrig,'Temperature','K');
title(sax(1),sprintf('Source window, Win Height: %.1f m, Win Dur: %.1f s, Win OL: %.2f'...
    ,winHeight*-D.dz,winDur*dt,tOver/winDur))

sax=plotThermStats(V0,V0,tTrig,'Velocity','m/s');
title(sax(1),sprintf('Source window, Win Height: %.1f m, Win Dur: %.1f s, Win OL: %.2f'...
    ,winHeight*-D.dz,winDur*dt,tOver/winDur))

% DETECTION WINDOW
dax=plotThermStats(Td,Vd,tTrig,'Temperature','K');
title(dax(1),sprintf('Detection window, Win Height: %.1f m, Win Dur: %.1f s, Win OL: %.2f'...
    ,round(winHeight.*trackParams.detection_window_height)*-D.dz,winDur*dt,tOver/winDur))

%% Plot all detection frames
% bpad = 15;
% imsz= size(squeeze(D.T(:,:,1)));
% 
% figure
% for pp=1:length(tI0)
%     tightSubplot(rows,cols,pp,0,0);
%     imagesc(squeeze(D.T(:,:,tI0(pp))))
%     
%     bbox0 = regionprops(D.mask(:,:,tI0(pp)),'BoundingBox');
%     bbox0 = bbox0.BoundingBox;
%     bbox = [max([0.5 bbox0(1)-bpad]) max([0.5 bbox0(2)-bpad]) bbox0(3:4)];
%     bbox(3:4) = [min([(imsz(2)+0.5)-bbox(1) bbox0(3)+bpad-(bbox(1)-bbox0(1))])...
%                   min([(imsz(1)+0.5)-bbox(2) bbox0(4)+bpad-(bbox(2)-bbox0(2))])];
%     axlims   = [bbox(1) sum(bbox([1 3])) bbox(2) sum(bbox([2 4])) ];
%     axis(axlims)
%     set(gca,'YDir','normal','XTick',[],'YTick',[])
%     caxis([240 350])
% end
% colormap(thermgray(150))

%% Check velocities
vI = 1;
nn = length(Vpos(vI).Vmu);
% Vmax = zeros(nn,1);
% Vmin = Vmax;
Vmax = 20;
Vmin = -10;
nbins = 60;
edges=linspace(Vmin,Vmax,nbins+1);
bins = edges(1:end-1)+diff(edges)/2;
hcounts = zeros(nn,nbins);
for fi=1:nn
    
    Vz = V.Vz(:,:,Vpos(vI).tI(fi));
    %     Vmax(fi) = max(Vz(Vpos(vI).pixelIdxList{fi}));
    %     Vmin(fi) = min(Vz(Vpos(vI).pixelIdxList{fi}));

    [hcounts(fi,:),~] = histcounts(Vz(Vpos(vI).pixelIdxList{fi}),edges);
    hcounts(fi,:) = hcounts(fi,:)./max(hcounts(fi,:));

%     histogram(Vz(Vpos(vI).pixelIdxList{fi}),edges);
%     title(sprintf('t = %.1f s',Vpos(vI).t(fi)))
%     pause(0.4)
end
figure
pcolor(Vpos(vI).t,bins,hcounts')
shading flat
hold on
plot(Vpos(vI).t,Vpos(vI).Vmu,'r','LineWidth',2)
%% Calc rise diagram
Rise = squeeze(max(D.T.*D.mask,[],2));
%% Plot rise diagram with tracks
figure
pcolor(D.t,D.z,Rise)
shading flat
% colormap(CubeHelix(200))
colormap(gray(200))
hold on
% Plot main windows
for ti = 1:length(Vpos)
    plotLineError(D.t(Vpos(ti).tI),[D.z(Vpos(ti).zI(1,:)) Vpos(ti).z D.z(Vpos(ti).zI(end,:))],[],0.4,true);
end
% Plot tracking windows
for ti = 1:length(Vpos)
    plotLineError(D.t(Vpos(ti).tI),[D.z(Vpos(ti).trackI(1,:)) D.z(Vpos(ti).trackI(end,:))],[0.85 0.325 0.1],0.2,true);
end
caxis([230 370])
xlabel('t [s]')
ylabel('z [m]')
ylim([40 450])

flargh
%% tracking/"vid" plots for a single source

% Plot window track w/ height
% figure
% pcolor(D.t,D.z,Tmax)
% shading flat
% colormap(thermgray(150))
% caxis([220 340])
% hold on
% % zlo = D.z( Vpos(vI0).zI(1,:) );
% % zhi = D.z( Vpos(vI0).zI(end,:) );
% % plotLineError(Vpos(vI0).t,[zhi Vpos(vI0).z zlo],[0 1 0.2],0.3,true)
% % for jj=1:numel(Vpos)
% %     zlo = D.z( Vpos(jj).zI(1,:) );
% %     zhi = D.z( Vpos(jj).zI(end,:) );
% %     plotLineError(Vpos(jj).t,[zhi Vpos(jj).z zlo],[0 1 0.2],0.3,true);
% % end
% set(gca,'FontSize',fs)
% xlim([0 175])
% ylim([min(D.z) 800])
% xlabel('Time [s]')
% ylabel('Height [m]')
% h=colorbar;
% h.Label.String = 'T_{max} [K]';

figure('position',[50 50 600 700])
for kk = 60 %[50 60 70] %1:numel(Vpos(vI0).t)
%     cla
    tI = Vpos(vI0).tI0+kk-1;
    zI = Vpos(vI0).zI(:,kk);
    z0 = T0.zI(:,1);
    srcadd = D.mask(:,:,tI);
    srcadd([1:z0(1)-1 z0(end)+1:end],:) = 0; % Zero outside window
    srcadd = D.T(:,:,tI).*edge(srcadd,'sobel');
 
    %     imgadd = D.T(:,:,tI).*D.mask(:,:,tI);
    imgadd = D.mask(:,:,tI);
    imgadd([1:zI(1)-1 zI(end)+1:end],:) = 0; % Zero outside window
    imgadd = D.T(:,:,tI).*edge(imgadd,'sobel');
%     srcadd = mask;
%     srcadd([1:S.zI(1,ii)-1 S.zI(end,ii)+1:end],:) = 0; % Zero outside window
    [cc,ch] = contour(D.x,D.z,srcadd,[0.9 0.9],'LineWidth',4,'Color',[0.850 0.325 0.098]);

   
    axes('position',[0.05 0.08 0.95 0.9])
    pcolor(D.x,D.z,D.T(:,:,tI)+srcadd) % +imgadd
    shading flat
    colormap(thermgray(150))
    caxis([230 340])
    daspect([1 1 1])
    axis([-300 300 40 800])
    set(gca,'FontSize',14)
    xlabel('x [m]')
    ylabel('z [m]')
    leg = sprintf('%i\n%.2f s',Vpos(vI0).tI0+kk-1,D.t(tI));
    t=text(0.8,0.8,leg,'FontSize',14,'Color', [0.9 .9 .9],'Units','normalized','HorizontalAlignment','right');
    
    fname = sprintf('feature_track_%i',Vpos(vI0).tI0+kk-1);
%     printpdf(fname,figdir,[15 18])
    
    pause(1)

end

%% Compare Source and tracked feature vs time
tStart = [D.t(Vpos(vI0).tI0+[50, 60, 70]-1)'; D.t(Vpos(vI0).tI0+[50, 60, 70]-1)'];
Trange1 = repmat([min(T0.prctile(:)); max(T0.prctile(:))], [1 size(tStart,2)] );
Trange2 = repmat([min(T(vI0).prctile(:)); max(T(vI0).prctile(:))], [1 size(tStart,2)] );


figure

% Temperatures
axa=tightSubplot(2,1,2);
plot(tStart(:,2),Trange1(:,2),'Color',[0.8 0.6 0.6],'LineWidth',2)
hold on
plotLineError(T0.t',T0.prctile');
% hold on
plot(T0.t,T0.mean,'Color', [0.85 0.325   0.098],'LineWidth',1.5)

legend('05-95','25-75','Median','Mean')
grid on
axis tight
set(gca,'FontSize',fs)
ylabel('T [K]')
xlabel('Time [s]')
% set(gca,'XTickLabel',[])
title(sprintf('Source history, Win Height: %.1f m, Win Dur: %.1f s, Win OL: %.2f'...
    ,winHeight*-D.dz,winDur*dt,tOver/winDur))

% Temperatures
axb=tightSubplot(2,1,1);
plot(tStart(:,2),Trange1(:,2),'Color',[0.8 0.6 0.6],'LineWidth',2)
hold on
plotLineError(T(vI0).t',T(vI0).prctile');
plot(T(vI0).t,T(vI0).mean,'Color', [0.85 0.325   0.098],'LineWidth',1.5)

% legend('05-95','25-75','Median','Mean')
grid on
axis tight
set(gca,'FontSize',fs)
ylabel('T [K]')
% title(sprintf('Tracked window, Win Height: %.1f m, Win Dur: %.1f s, Win OL: %.2f'...
%     ,winHeight*-D.dz,winDur*dt,tOver/winDur))
text(0.75, 0.8, 'Tracked Window','FontSize',14,'Units','normalized')
% xlabel('Time [s]')
set(gca,'XTickLabel',[])

linkaxes([axa axb],'xy')

%% Quick and dirty exponential
% Initial mixture
phi_v = 0.1;       % Gas mass fraction
phi_m = 1-phi_v;    % Particle mass fraction
rho_v = 0.2;        % Gas density
rho_m = 2500;       % Particle density

[ep,rho_p] = wt2vol([phi_m phi_v],[rho_m rho_v]);

C_m = 1300; % Heat capacity, ash
C_v = 1996; % Heat capacity, vapour
C_p0 = phi_m*C_m + phi_v*C_v;

C_a = 1004; % Heat capacity, air

% Overwrite
% C_p0 = C_a;
rho_p = 1;

R      =  [10 100];
Vrise   = 10; 
alpha   = [0.1 0.25];
C_p     = 1000;
rho_a   = 1;
C_a     = 1004;


R0 = 10;
m0 = 4/3*pi*R0.^3*rho_p;
Gamma = m0*C_p; % Constant Gamma (wrong, but hey...)
Ve      = Vrise*alpha(1);

iStart = 1;
dT0=sqrt(T(vI0).var(iStart));
n = 1000;
tf = 50;
t1 = linspace(0,tf,n); % Seconds, for now

dT1 = dT0*exp(-Ve.^3./3*(rho_p*C_p/Gamma).*t1.^3);
efold = (1/Ve.^3./3*(rho_p*C_p/Gamma)).^(1/3);
% 
figure

plot(t1,dT1/dT0);
xlabel('t [s]')
title('Fixed CV solution')



%% Plot T stats for tracked feature(s) - vs time

t = Vpos(vI0).t-Vpos(vI0).t(1);


figure
% Velocity distributions
axa=tightSubplot(3,1,1,[],0.01);
rectangle('Position',[t(1) 0 19 15],'FaceColor',[0.95 0.9 0.8 0.3],'EdgeColor','none')
hold on
for ii=1:6; plot(Vpos(ii).t-Vpos(ii).t(1),smooth(Vpos(ii).Vmu,20),'LineWidth',0.8,'Color',[0.4 0.4 0.4]); hold on; end
% plot(t,Vpos(vI0).Vmu,'Color', rgba2rgb([0   0.447   0.741],0.7,[1 1 1]));
axis tight
% hold on
plot(t,smooth(Vpos(vI0).Vmu,20),'Color', [0   0.447   0.741],'LineWidth',2)
grid on
set(gca,'FontSize',fs)
set(gca,'XTickLabel',[])

% pcolor(D.t,Tbins,HistC0); shading flat
ylabel('v_y [m/s]')
% title(sprintf('Tracked window, Win Height: %.1f m, Win Dur: %.1f s, Win OL: %.2f'...
%     ,winHeight*-D.dz,winDur*dt,tOver/winDur))

% Temperatures
axb=tightSubplot(3,1,2,[],0.01);
% rectangle('Position',[11 250 19 142],'FaceColor',[0.95 0.9 0.8 0.3],'EdgeColor','none')
% hold on
plotLineError(t,T(vI0).prctile');
hold on
plot(t,T(vI0).mean,'Color', [0.85 0.325   0.098],'LineWidth',1.5)
legend('05-95','25-75','Median','Mean')
grid on
axis tight
set(gca,'FontSize',fs)
ylabel('T_b [K]')
set(gca,'XTickLabel',[])
ylim([250 410])

% Temperature variance
% axb=tightSubplot(3,1,2);
% plot(T(vI0).t,T(vI0).var,'LineWidth',1.5)
% ylabel('T variance [K^2]')
% grid on
% axis tight
% set(gca,'FontSize',12)
% set(gca,'XTickLabel',[])

% Temperature variance
axc=tightSubplot(3,1,3,[],0.01);
% rectangle('Position',[11 0 17.5 15],'FaceColor',[0.8 0.8 0.8 0.3],'EdgeColor','none')
% hold on
plot(t,sqrt(T(vI0).var),'LineWidth',1.5)
hold on
% plot(t1+T(vI0).t(1),dT1);
ylabel('T_b Std. Dev. [K]')
grid on
axis tight
set(gca,'FontSize',fs)



xlabel('Time [s]')
linkaxes([axa axb axc],'x')
xlim([0.5 40])

%% Plot...something about the tracks?
figure
plotLineError(t,T(vI0).prctile');
hold on
pa=plot(t,T(vI0).mean,'Color', 'k','LineWidth',2);
pb=plot(t,T(vI0).prctile(5,:), 'k--','LineWidth',2);
legend([pa pb],{'Mean','95th percentile'})
grid on
axis tight
set(gca,'FontSize',fs)
ylabel('T_b [K]')
xlabel('t [s]')
% set(gca,'XTickLabel',[])
ylim([250 410])
xlim([0 25])
%% Plot T stats for tracked feature(s) - vs height
figure
atmoT = interp1(atmoZ,atmo.Temperature,T(vI0).z,'linear')-20;

% Temperatures
axa=tightSubplot(3,1,1);
plotLineError(T(vI0).z',T(vI0).prctile');
% plotLineError(T(vI0).z',T(vI0).prctile'-atmoT');
hold on
plot(T(vI0).z,T(vI0).mean,'Color', [0.85 0.325   0.098],'LineWidth',1.5)
% plot(atmo.Height-ventz,atmo.Temperature,'--k')
legend('05-95','25-75','Median','Mean')
grid on
axis tight
set(gca,'FontSize',12)
ylabel('T [K]')
set(gca,'XTickLabel',[])
title(sprintf('Tracked window, Win Height: %.1f m, Win Dur: %.1f s, Win OL: %.2f'...
    ,winHeight*-D.dz,winDur*dt,tOver/winDur))

% Temperature variance
axb=tightSubplot(3,1,2);
plot(T(vI0).z,T(vI0).var,'LineWidth',1.5)
ylabel('T variance [K^2]')
grid on
axis tight
set(gca,'FontSize',12)
set(gca,'XTickLabel',[])

% Velocity distributions
axc=tightSubplot(3,1,3);
plot(Vpos(vI0).z,Vpos(vI0).Vmu);
axis tight
hold on
plot(Vpos(vI0).z,smooth(Vpos(vI0).Vmu,20),'Color', [0.85 0.325   0.098],'LineWidth',1.5)
grid on
set(gca,'FontSize',12)
% pcolor(D.t,Tbins,HistC0); shading flat
ylabel('v_y [m/s]')
xlabel('Height [m]')
linkaxes([axa axb axc],'x')

% xlim([min(T(vI0).z) max(T(vI0).z)])

%% Variances for all tracks
% Get some non-dimensionalized times and heights?
tmax = 100;
tlimsV=[3 12];
vmu = zeros(length(Vpos));
for kk=1:length(Vpos)
    tv  = Vpos(vI0).t-Vpos(vI0).t(1);
    tvI = tv<100;
    vmu(kk) = mean(Vpos(kk).Vmu(and(tv<tlimsV(2),tv>tlimsV(1))));
end

L = 75; %100?  % Characteristic length (plume diameter?)
tprime = L./vmu;
tprimeLims = [5 9];

%%
R      =  [25];
alpha  =  [0.1 0.25];
Tstart =  1; %T(vI0).prctile(1,5);
[T2,t2,Vrise] = mixingModel(Tstart,100,alpha,R);



figure
plot(t2/(L/Vrise)-.6,T2(:,1),'-.','Color',[0.5 0.5 0.5],'LineWidth',2);
hold on
plot(t2/(L/Vrise)-.4,T2(:,2),'--','Color',[0.75 0.75 0.75],'LineWidth',2);

for kk=1:numel(T)

%     SigSmooth = smooth(sqrt(T(kk).var./Tm),20);
    Sig = sqrt(T(kk).prctile(5,:));
% tightSubplot(2,1,1);
%     [Tm,ii] = max(T(kk).var);
    [Sm,ii] = max(Sig);
    tSigSmooth = (T(kk).t-T(kk).t(ii))/tprime(kk);
    
    SigBase = mean(Sig(and(tSigSmooth>tprimeLims(1),tSigSmooth<tprimeLims(2))));
    SigSmooth = smooth((Sig-SigBase)/(Sm-SigBase),20);
    
    plot(tSigSmooth,SigSmooth,'LineWidth',1.5)
%     plot((T(kk).t-T(kk).t(ii))/tprime(kk),smooth(sqrt(T(kk).var./Tm),20),'LineWidth',1.5)
%     plot((T(kk).t-T(kk).t(ii)),sqrt(T(kk).var./Tm),'LineWidth',1.5)
    hold on
end
ylabel('$\frac{\Delta T}{\Delta T_0}$','interpreter','latex')
xlabel('$t/t_{eddy}$','interpreter','latex')

hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right')
grid on
axis tight
set(gca,'FontSize',14)
xlim([0 8])
% xlabel('Height [m]')
legend('alpha=0.1','alpha=0.25','Pulse 1','Pulse 2','Pulse 3','Pulse 4','Pulse 5')
% 
% tightSubplot(2,1,2);
% for kk=1:numel(T)
% %     [Tm,ii] = max(T(kk).var);
%     plot(Vpos(kk).z,smooth(Vpos(kk).Vmu,20),'LineWidth',1.5)
%     hold on
% end
% ylabel('v [m/s]')
% xlabel('Height [m]')
% xlim([0 400])

%% Lil bit o slice
% Ts = D.T(1:250,1:200,1:2:600);
% Ts = permute(fliplr(Ts), [2 3 1]);
% 
% figure
% % slice(D.x(1:300),D.z(1:250),D.t(1:2:600),Ts,50,50,50)
% slice(D.t(1:2:600),D.x(1:200),D.z(1:250),Ts,100,150,50)
% shading flat
% xlabel('t [s]')
% ylabel('x [m]')
% zlabel('z [m]')
% daspect([0.1 1 1])
% colormap(thermgray(150))
% caxis([230 340])
% set(gca,'FontSize',fs)

% Time vertical
% hslice = squeeze(D.T(10,:,:));
% 
% figure
% pcolor(D.x(1:200),D.t(1:700),hslice(1:200,1:700)')
% shading flat
% xlabel('x [m]')
% ylabel('t [s]')
% zlabel('z [m]')
% daspect([1 0.15 1])
% colormap(thermgray(150))
% caxis([230 340])
% set(gca,'FontSize',fs)
% view([0 0 1])
% axis tight


%% Velocity image
fi = 300;
figure('position',[50 50 600 700])

pcolor(D.x,D.z,D.T(:,:,fi)+D.T(:,:,fi).*edge(D.mask(:,:,fi),'sobel')); shading flat
colormap(thermgray(150))
caxis([230 340])
hold on
% plotThermVelocities(D.x,D.z,V.Vx.*D.mask,V.Vz.*D.mask,15,fi,'vector') %,pt,ROI)
xlim([-315 400])
ylim([50 900])
set(gca,'FontSize',12)
xlabel('x [m]')
ylabel('z [m]')
daspect([1 1 1])
% axis equal
% cbar = colorbar(


%% Source conditions with horizontal cube slice
% Time Horizontal
hslice = squeeze(D.T(10,:,:));

figure('position',[50 50 900 700])
% set(gcf,'Color','none')
% set(gcf,'InvertHardCopy','off')

axa=tightSubplot(2,1,1,[],[],[],[],[1.5 1 1]);
pcolor(D.t(1:700),D.x(1:200),hslice(1:200,1:700))
shading flat
ylabel('x [m]')
xlabel('t [s]')
zlabel('z [m]')
% daspect([0.15 1 1])
colormap(thermgray(150))
caxis([230 340])
set(gca,'FontSize',fs)
set(gca,'XTickLabel',[])
view([0 0 1])
axis tight
ylim([-225 75])
xlim([0 150])
set(gca,'XColor',fc,'YColor',fc)
% tightSubplot(3,1,2)
cb = colorbar('north');
cb.FontSize=14;
% cb.Color=[0.85 0.85 0.85];
cb.Position = [0.62 0.87 0.3 0.05];
cb.Label.String='T_b [K]';
cb.Label.Position=[220 1];

% c.Label.Position=[0.6 0.9];
cb.AxisLocation='out';
% cb.Label.Color=[1 1 1]; 
% grid on
% set(gca,'GridColor',[0.8 0.8 0.8])

% Temperatures
axb=tightSubplot(2,1,2,[],[],[],[],[1.5 1 1]);
plotLineError(T0.t',T0.prctile');
hold on
plot(T0.t,T0.mean,'Color', [0.85 0.325   0.098],'LineWidth',1.5)
legend('05-95','25-75','Median','Mean')
grid on
axis tight
set(gca,'FontSize',fs)
% set(gca,'XTickLabel',[])
ylabel('T [K]')
xlim([0 150])
set(gca,'XColor',fc,'YColor',fc)

% set(gca,'XTickLabel',[])
% title(sprintf('Source history, Win Height: %.1f m, Win Dur: %.1f s, Win OL: %.2f'...
%     ,winHeight*-D.dz,winDur*dt,tOver/winDur))

% Temperature variance
% axc=tightSubplot(3,1,3,[],[],[],[],[1.5 1 1]);
% plot(T0.t,T0.var,'LineWidth',1.5)
% hold on
% % plot(tStart,Trange,'Color',[0.8 0.6 0.6],'LineWidth',2.5)
% 
% ylabel('T variance [K^2]')
% grid on
% axis tight
xlabel('Time [s]')
% set(gca,'FontSize',fs)
% xlim(axc,[0 150])
% grid on

% set(gca,'XTickLabel',[])

% Velocity distributions
% axc=tightSubplot(3,1,3);
% plotLineError(V0.t',V0.prctile');
% axis tight
% hold on
% plot(V0.t,smooth(V0.mean,10),'Color', [0.85 0.325   0.098],'LineWidth',1.5)
% % grid on
% set(gca,'FontSize',fs)
% % pcolor(D.t,Tbins,HistC0); shading flat
% ylabel('v_y [m/s]')
% xlabel('Time [s]')
% linkaxes([axa axb axc],'x')
% linkaxes([axa axb axc],'x')

%% Plot mean image and profiles against pulse profiles
% Mean image
fi = 300;
figure('position',[50 50 650 700])

axes('position',[0.05 0.08 0.94 0.88])
pcolor(Dmu.x,Dmu.z,Dmu.T+Dmu.T.*edge(Dmu.mask,'sobel')); shading flat
colormap(thermgray(150))
caxis([min(Dmu.T(:)) max(Dmu.T(:))])
hold on
xlim([-315 400])
ylim([50 900])
set(gca,'FontSize',14)
xlabel('x [m]')
ylabel('z [m]')
daspect([1 1 1])
cb = colorbar('west');
cb.Position= [0.1467    0.65    0.0333    0.3];
cb.FontSize=14;
cb.Color=[0.92 0.92 0.92];
cb.Label.String='T_b [K]';
cb.Label.Color=[0.92 0.92 0.92];  
cb.AxisLocation='in';

% Compare height profiles of blobs vs mean image
figure
plot(Tmu.z,Tmu.prctile(5,:),'Color',[0   0.447   0.741],'LineWidth',2)
hold on
% plot(Tmu.z,Tmu.prctile(3,:),'Color',rgba2rgb([0   0.447   0.741],0.5,get(gca,'Color')),'LineWidth',2)
% plot(Tmu.z,Tmu.mean,'Color',rgba2rgb([0   0.447   0.741],0.5,get(gca,'Color')),'LineWidth',2)

for ii=1:3; %length(T)
    plot(T(ii).z,T(ii).prctile(5,:),'Color',[0.850 0.325 0.098],'LineWidth',1)
%     plot(T(ii).z,T(ii).prctile(3,:),'Color',rgba2rgb([0.850 0.325 0.098],0.5,get(gca,'Color')),'LineWidth',1)
%     plot(T(ii).z,T(ii).mean,'Color',rgba2rgb([0.850 0.325 0.098],0.5,get(gca,'Color')),'LineWidth',1)
end

xlim([0 400])
xlabel('z [m]')
xlabel('T_b [K]')

% Normalized height profiles of blobs vs mean image
Tmu_av = Tmu.prctile(5,:);
Tmu_z = Tmu.z(Tmu.z<=400);
Tmu_av = Tmu_av(Tmu.z<=400);

figure
plot(Tmu_z,(Tmu_av-min(Tmu_av))/range(Tmu_av),'Color',[0   0.447   0.741],'LineWidth',2)
hold on
% plot(Tmu.z,Tmu.prctile(3,:),'Color',rgba2rgb([0   0.447   0.741],0.5,get(gca,'Color')),'LineWidth',2)
% plot(Tmu.z,Tmu.mean,'Color',rgba2rgb([0   0.447   0.741],0.5,get(gca,'Color')),'LineWidth',2)

for ii=1:length(T)
    T_blob = T(ii).prctile(5,:);
    T_z = T(ii).z(T(ii).z<=400);
    T_blob = T_blob(T(ii).z<=400);
    plot(T_z,(T_blob-min(T_blob))/range(T_blob),'Color',[0.850 0.325 0.098],'LineWidth',1)
%     plot(T(ii).z,T(ii).prctile(3,:),'Color',rgba2rgb([0.850 0.325 0.098],0.5,get(gca,'Color')),'LineWidth',1)
%     plot(T(ii).z,T(ii).mean,'Color',rgba2rgb([0.850 0.325 0.098],0.5,get(gca,'Color')),'LineWidth',1)
end

xlim([min(Tmu.z) 400])
xlabel('z [m]')
ylabel('T_b [K]')

%% All tracked velocities
figure('Position',[200 100 800 350])
rectangle('Position',[125 0 200 15],'FaceColor',[0.8 0.8 0.8 0.3],'EdgeColor','none')
hold on
for ii=1:6; plot(Vpos(ii).z,smooth(Vpos(ii).Vmu,20),'LineWidth',1.5); hold on; end
grid on
xlim([100 475])
set(gca,'FontSize',14)
ylabel('u [m/s]')
xlabel('z [m]')
% -Vpos(ii).z(1)
%% Spit out simplified tracks
 for ii=1:length(Vpos)
     Ttracks(ii).t=Vpos(ii).t; 
     Ttracks(ii).z=Vpos(ii).z; 
     Ttracks(ii).Vavg=smooth(Vpos(ii).Vmu,20); 
     Ttracks(ii).Tavg=T(ii).mean'; 
     Ttracks(ii).Tvariance=T(ii).var'; 
     Ttracks(ii).T95 = T(ii).prctile(5,:)'; 
     
     Track=struct2table(Ttracks(ii));
     writetable(Track,fullfile(cubeDir,sprintf('Track%i.csv',ii)))
%      save(fullfile(cubeDir,sprintf('Track%i',ii)),'Track','-asii','-tabs')
%      Track
 end
save(fullfile(cubeDir,'Ttracks'),'Ttracks')

%% Functions but something may have been deleted...
% function [idx,t] = getSTFTColumns(nx,nwin,noverlap,Fs)
% % Borrowed from the MATLAB spectrogram function
% % IN:
% % x        = input signal
% % nx       = length of input signal
% % nwin     = length of each window
% % noverlap = numner of samples each segment overlaps
% % Fs       = sampling freq
% %
% % OUT:
% % idx      = array indices
% % Determine the number of columns of the STFT output (i.e., the S output),
% % the times, t centered on windows, and the associated matrix indices
% ncol = fix((nx-noverlap)/(nwin-noverlap));
% 
% colindex = 1 + (0:(ncol-1))*(nwin-noverlap);
% rowindex = (1:nwin)';
% % 'xin' should be of the same datatype as 'x'
% % xin = zeros(nwin,ncol,class(x)); %#ok<*ZEROLIKE>
% 
% % Put x into columns of xin with the proper offset
% idx = rowindex(:,ones(1,ncol))+colindex(ones(nwin,1),:)-1;
% % xin(:) = x(idx);
% 
% % colindex already takes into account the noverlap factor; Return a T
% % vector whose elements are centered in the segment.
% t = ((colindex-1)+((nwin)/2)')/Fs;
% end