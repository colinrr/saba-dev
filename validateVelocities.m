% =====================================================================
%          Comparing image velocities
% =====================================================================
% Run test calibration of optical flow velocities against those retrieved
% from plume masks and from cross-correlation

 clear all; close all;
%% USER INPUT

% -----------DATA DIRECTORIES------------
dataDir = '~/Kahuna/data/sabancaya_5_2018/image_exports/24A/';
interpDir = fullfile(dataDir,'interp-mat-test/');

dataCube  = fullfile(dataDir, 'thermCubeAnalysis/thermStats_2019-09-18_z641_x591_t1195.mat');
dataCube2 = fullfile(dataDir, 'thermCubeAnalysis/thermStats_2019-07-05_z131_x161_t598.mat');


ovid = fullfile(dataDir,'vids/thermal-gradient'); % No extension

velCube = fullfile(dataDir, 'thermCubeAnalysis/velocimetry_190918_n1194_nSz5_nPyr3_fSz15.mat');

thermCorr = fullfile(dataDir, 'thermCubeAnalysis/thermCorrelation_2019-07-05_z131_x161_t598.mat');
% Compare 3 velocities: initial plume top, xcorr time lag, and opticFlow

%% (1) Get mask velocity
load(dataCube)

N = numel(D.idx);
maskLine = squeeze(sum(D.mask,2));
t = D.t;
%% 
zMask = zeros(N,1);
xMask = zMask;
hJ = ones(N,1);
hI = hJ;

for ii=1:N
    htopJ = find(maskLine(:,ii)>0,1,'last');
    if ~isempty(htopJ)
        zMask(ii) = D.z(htopJ);
        hJ(ii) = htopJ;
        
        htopI = round( mean(find(D.mask(htopJ,:,ii))));
        if ~isempty(htopI)
            hI(ii) = htopI;
        end
    end
end

zMask0 = D.z(1);
zMask = zMask-zMask0;
zMask = smooth(zMask,5);
spoints = 20;
[vfM,vsM,tv] = dHdt(D.t,zMask,spoints);

% clear D
%% (2) Get velocities using mask positioning from above and OpticFlow values
load(velCube)

%%
rD = 11; % RMS vels of a region this wide/high in pixels

vxOF = zeros(N-1,1);
vzOF = zeros(N-1,1);
vxOFr = zeros(N-1,1);
vzOFr = zeros(N-1,1);

for ii=1:N-1
    % Using just a single value
    vxOF(ii) = V.Vx(hJ(ii), hI(ii), ii); 
    vzOF(ii) = V.Vz(hJ(ii), hI(ii), ii);
    
    % Consider using an RMS value of local V's (5x5?)
    hJrms = hJ(ii)-(rD-1):hJ(ii);
    hIrms = hI(ii)-floor(rD/2):hI(ii)+floor(rD/2);
    
    if any(hJrms(:)<1)
        hJrms = hJrms - min(hJrms(:))+1;
    elseif any(hJrms(:)>size(V.Vx,1))
        hJrms = hJrms - (max(hJrms)-size(V.Vx,1));
    end    
    
    if any(hIrms(:)<1)
        hIrms = hIrms - min(hIrms(:))+1;
    elseif any(hIrms(:)>size(V.Vx,2))
        hIrms = hIrms - (max(hIrms)-size(V.Vx,2));
    end
    
    
    vxt = V.Vx(hJrms,hIrms,ii);
    vzt = V.Vz(hJrms,hIrms,ii);
    vxOFr(ii) = rms(vxt(:)); 
    vzOFr(ii) = rms(vzt(:));
    
 
end

vzOFs  = smooth(vzOF,20);
%% Velocities from xcorr lags
load(thermCorr)
load(dataCube2)
dz = mean(diff(D.z));
dt = mean(diff(D.t));

zCorr = D.z;
tCorr = -dt*C(1).Lag0;
tCorr = tCorr - min(tCorr);

% t2 = linspace(0, max(tCorr),numel(tCorr));


tCorr = smooth(tCorr,5);
spoints = 5;
[vfC,vsC,tvC] = dHdt(tCorr,zCorr,spoints);

%% Plot em up
figure
ax1=tightSubplot(3,1,1,[],0);
plot(t,zMask,'.-','LineWidth',2)
ylabel('H [m]')
title('Plume-top velocity comparisons')
set(gca,'FontSize',12)

ax2=tightSubplot(3,1,2,[],0);
plot(tv,vfM)
hold on
plot(t(2:end),vzOFr) % Vertical component only
plot(tvC,vfC)
ylabel('Raw v [m/s]')
set(gca,'FontSize',12)


ax3=tightSubplot(3,1,3,[],0);
plot(tv,vsM,'LineWidth',2)
hold on
plot(t(2:end),vzOFs,'LineWidth',2)
plot(tvC,vsC,'LineWidth',2)
ylabel('Smoothed v [m/s]')
legend('Mask top','OpticFlow','XCorr')
linkaxes([ax1 ax2 ax3],'x')
axis tight
set(gca,'FontSize',12)

xlabel('t [s]')

%%
figure
plot(tv,vsM,'LineWidth',2)
hold on
plot(t(2:end),vzOFs,'LineWidth',2)
plot(tvC,vsC,'LineWidth',2)
ylabel('U [m/s]')
legend('Mask top','OpticFlow','XCorr')
linkaxes([ax1 ax2 ax3],'x')
axis tight
set(gca,'FontSize',14)

xlabel('t [s]')
xlim([0.5 100])
%% Velocity function
function [vf,vs,tv] = dHdt(t,h,s)
% Estimates velocities from time,height data using FD scheme
% s  = length of moving average, defaults to 5
%
% vf = (non-uniform?) finite difference scheme
% vs = smoothed with an s-point average
% tv = output times

if nargin<3
    s = 5;
end

% Plume tracker start with the second frame, so let's assume a t and h
% start for now
% t = [0;t]; % Assumed point only accurate to first time step
% h = [0;h];

% For now, let's quickly interpolate to regular time stamps
% Later, we'll cut out anything with timestamps too large
 %  - choosing just 1 second for now
dt = diff(t);
if any(dt<1)
    dt0 = round(mean(dt(dt<1)),2);
else
    dt0 = round(mean(dt));
end
tn = [t(1):dt0:t(end)]';
hn = interp1(t,h,tn,'linear');
% hn = smooth(hn,5);

n = length(hn); 
e = ones(n-1,1);
% Derivative matrix, nodes to cells
Dn2c = 1/dt0 * spdiags([-e,e],[0,1],n-1,n); 
% Averaging matrix, cells to nodes, ignoring boundary values atm
Ac2n = 1/2   * spdiags([e,e],[0,1],n-2,n-1); 

% dt = diff(t);
% t_d = dt/2+t(1:end-1);

% dti = dt(2:end); dti_1 = dt(1:end-1);
% dt_1 = -dti./(dti_1.*(dti-dti_1));
% dt_2 = (dti-dti_1)./(dti.*dti_1);
% dt_3 = dti_1./(dti.*(dti+dti_1));
% D = spdiags([dt_1 dt_2 dt_3],[0,1,2],n-2,n-2);

% Run the finite difference
vf = Ac2n * Dn2c * hn;
vs = smooth(vf,s);
tf = tn(2:end-1);

% Interpolate back to frame times?
tv = t; % Match time vector for now
vf = interp1(tf,vf,tv,'linear',NaN);
vs = interp1(tf,vs,tv,'linear',NaN);


% vl = diff(h)./dt; % Results in a cell-centered deal
% vf = D*h(2:end-1);

end

%% Deprecated junkpile...
% 
% opticIdx = []; %[400:2:548];
% % -----------VIDEO PARAMS------------
% 
% vidPar.ROI      = []; %[370 600 580 755];
% vidPar.Tthresh  = 240;
% vidPar.Gthresh  = 150; % Maximum temperature gradient - values will be normalized to this
%                        % At the moment, only able to manually pick this
%                        % from the hottest image
% vidPar.gamma    = 0.5;
% vidPar.fgmaskFile = fullfile(dataDir,'interp-mat/foreground_mask_edit.mat');
% 
% % Vid flags
% plumeMask  = false;
% scaleFrame = true;
% plotFrames = true;
% 
% vidFlags = [plumeMask scaleFrame plotFrames];
% % -----------OPTICAL FLOW------------
% % vidParFile = fullfile(dataDir,'vids/thermal-gradient_params2019-09-08_75frames.mat');
% 
% geomf = fullfile(dataDir,'mat/geometry.mat');
% T = load(geomf,'T');
% T = T.T;
% 
% % VelFile = 
% % -----------DRIVER FLAGS------------
% flag_makeVid   = false;
% flag_opticFlow = false;
% 
% % DO THE THINGS
% 
% if flag_makeVid
%     [vidParams,vidParFile] = thermal2gradient(interpDir,T,opticIdx,vidPar,vidFlags,ovid);
% end
% 
% if flag_opticFlow
% %     if ~isempty(vidParFile)
% %         load(vidParFile)
% %     end
%     [Vx,Vz] = thermVelocity(ovid, vidParFile);
% end
