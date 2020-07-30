% Compare two velocity sets - original thermOpticFlow vs thermOpticFlow2O
clear all; close all

homeDir = '~'; % LINUX/MAC
dataDir   = fullfile(homeDir,'/Kahuna/data/sabancaya_5_2018/');

% -------- MAIN WORKING directory for this event --------
thermDir   = fullfile(dataDir,'image_exports/25B/');
cubeDir = fullfile(thermDir,'thermCubeAnalysis/');

thermFile   = fullfile(cubeDir,'thermStats_2020-03-25_z710_x571_t1294.mat');
velCube1    = fullfile(cubeDir, 'opticFlowCNL_20-03-26_n1293.mat');
velCube2    = fullfile(cubeDir, 'opticFlowCNL_20-04-09_n1294.mat');

D  = loadif(thermFile,'D');
V1 = loadif(velCube1,'V');
V2 = loadif(velCube2,'V');


%% Show all-plume velocity/angle average/std changes
Vmu1 = zeros(length(V1.t),1);
Vsig1 = Vmu1;
Umu1 = Vmu1;
Usig1 = Vmu1;
for ii=19:length(Vmu1)-1
    Vx = V1.Vx(:,:,ii);
    Vz = V1.Vz(:,:,ii);
    
    Uabs = sqrt(Vx.^2 + Vz.^2);
    Theta = atan2d(Vz,Vx);
    
    % Average and standard deviation in U,V for each plume mask
    Vmu1(ii) = mean(Vz(D.mask(:,:,ii)));
    Vsig1(ii) = std(Vz(D.mask(:,:,ii)));
    Umu1(ii) = mean(Vx(D.mask(:,:,ii)));
    Usig1(ii) = std(Vx(D.mask(:,:,ii)));
    
    % Average and standard deviation in abs(U) and vector angle Theta for
    % each plume mask
    UabsMu1(ii) = mean(Uabs(D.mask(:,:,ii)));
    UabsSig1(ii) = std(Uabs(D.mask(:,:,ii)));
    ThetaMu1(ii) = mean(Theta(D.mask(:,:,ii)));  
    ThetaSig1(ii) = std(Theta(D.mask(:,:,ii)));

end


Vmu2 = zeros(length(V2.t),1);
Vsig2 = Vmu2;
Umu2 = Vmu2;
Usig2 = Vmu2;
for ii=19:length(Vmu1)
    Vx = V2.Vx(:,:,ii);
    Vz = V2.Vz(:,:,ii);
    
    Uabs = sqrt(Vx.^2 + Vz.^2);
    Theta = atan2d(Vz,Vx);
    
    Vmu2(ii) = mean(Vz(D.mask(:,:,ii)));
    Vsig2(ii) = std(Vz(D.mask(:,:,ii)));
    Umu2(ii) = mean(Vx(D.mask(:,:,ii)));
    Usig2(ii) = std(Vx(D.mask(:,:,ii)));
    
    UabsMu2(ii) = mean(Uabs(D.mask(:,:,ii)));
    UabsSig2(ii) = std(Uabs(D.mask(:,:,ii)));
    ThetaMu2(ii) = mean(Theta(D.mask(:,:,ii)));  
    ThetaSig2(ii) = std(Theta(D.mask(:,:,ii)));

end

%%  Calc and plot total averages/std for Vz 
Vmumu1 = nanmean(Vmu1);
Vsigmu1 = nanmean(Vsig1);
Vmusig1 = nanstd(Vmu1);

Vmumu2 = nanmean(Vmu2);
Vsigmu2 = nanmean(Vsig2);
Vmusig2 = nanstd(Vmu2);

[Vmumu1 Vmusig1; Vmumu2 Vmusig2]

figure
tightSubplot(2,1,1)
% plotLineError(V1.t(1:end-1),[Tmumu1-Tsigmu1 Tmumu1 Tmumu1+Tsigmu1])
plot(V1.t(1:end-1)+diff(V1.t)/2,Vmu1(1:end-1))
hold on
plot(V2.t,Vmu2)
tightSubplot(2,1,2)
% plotLineError(V2.t,[Tmumu2-Tsigmu2 Tmumu2 Tmumu2+Tsigmu2])
plot(V1.t(1:end-1)+diff(V1.t)/2,Vsig1(1:end-1))
hold on
plot(V2.t,Vsig2)

%% % Calc average/std for absolute velocity and direction
Uabsmumu1 = nanmean(UabsMu1);
Uabsmusig1 = nanstd(UabsMu1);
ThetaMumu1 = nanmean(ThetaMu1);
ThetaMusig1 = nanstd(ThetaMu1);

Uabsmumu2 = nanmean(UabsMu2);
Uabsmusig2 = nanstd(UabsMu2);
ThetaMumu2 = nanmean(ThetaMu2);
ThetaMusig2 = nanstd(ThetaMu2);

figure
tightSubplot(2,1,1)
% plotLineError(V1.t(1:end-1),[Tmumu1-Tsigmu1 Tmumu1 Tmumu1+Tsigmu1])
plot(V1.t(1:end-1)+diff(V1.t)/2,UabsMu1)
hold on
plot(V2.t,UabsMu2)
tightSubplot(2,1,2)
% plotLineError(V2.t,[Tmumu2-Tsigmu2 Tmumu2 Tmumu2+Tsigmu2])
plot(V1.t(1:end-1)+diff(V1.t)/2,UabsSig1)
hold on
plot(V2.t,UabsSig2)


figure
tightSubplot(2,1,1)
% plotLineError(V1.t(1:end-1),[Tmumu1-Tsigmu1 Tmumu1 Tmumu1+Tsigmu1])
plot(V1.t(1:end-1)+diff(V1.t)/2,ThetaMu1)
hold on
plot(V2.t,ThetaMu2)
tightSubplot(2,1,2)
% plotLineError(V2.t,[Tmumu2-Tsigmu2 Tmumu2 Tmumu2+Tsigmu2])
plot(V1.t(1:end-1)+diff(V1.t)/2,ThetaSig1)
hold on
plot(V2.t,ThetaSig2)

%% Abs velocity and direction a few single pixels

% iPix = [100 160 260 75];
% jPix = [190 100 100 290];

iPix = [100 75];
jPix = [190 290];

II = sub2ind([length(V1.z) length(V1.x)],iPix,jPix);

U1abs_t = zeros(size(V1.Vz,3),length(iPix));
Theta1_t = U1abs_t;
t1 = V1.t(1:end-1)+diff(V1.t)/2;

U2abs_t = zeros(size(V2.Vz,3),length(iPix));
Theta2_t = U1abs_t;

for ii=19:size(U1abs_t,1)
    Vx = V1.Vx(:,:,ii);
    Vz = V1.Vz(:,:,ii);
    
    U1abs_t(ii,:) = sqrt(Vx(II).^2 + Vz(II).^2);
    Theta1_t(ii,:) = atan2d(Vz(II),Vx(II));
    
    U1std = movstd(U1abs_t,100,1);
end
for ii=19:size(U2abs_t,1)
    Vx = V2.Vx(:,:,ii);
    Vz = V2.Vz(:,:,ii);
    
    U2abs_t(ii,:) = sqrt(Vx(II).^2 + Vz(II).^2);
    Theta2_t(ii,:) = atan2d(Vz(II),Vx(II));
    
    U2std = movstd(U2abs_t,100,1);
end

co = get(0,'DefaultAxesColorOrder');
co0 = rgb2hsv(co);
co0(:,2) = co0(:,2)*0.25;
co0 = hsv2rgb(co0);

% Plot velocity and angle vectors
figure
axb(1)=tightSubplot(2,1,1);
% ax(1).ColorOrder = [co0(1:length(iPix),:) ; co(1:length(iPix),:)];
hold on
for ii=1:length(iPix)
    plot(t1,U1abs_t(:,ii),'LineWidth',1.5,'Color',co0(ii,:));
end
for ii=1:length(iPix)
    plot(V2.t,U2abs_t(:,ii),'LineWidth',1.2,'Color',co(ii,:));
end
legend('Inside mask, 1st Order','Outside mask, 1st Order','Inside mask, 2nd Order','Outside mask, 2nd Order')
title('Values of two pixels over time, 1st- vs 2nd-order Optical Flow')
set(gca,'FontSize',12)
ylabel('|U| [m/s]')

axb(2)=tightSubplot(2,1,2);
for ii=1:length(iPix)
    plot(t1,Theta1_t(:,ii),'LineWidth',1.5,'Color',co0(ii,:));
    hold on
end
for ii=1:length(iPix)
    plot(V2.t,Theta2_t(:,ii),'LineWidth',1.2,'Color',co(ii,:));
end
xlabel('Time [s]')
set(gca,'FontSize',12)
ylim([-180 180])
ylabel('Vector Angle [degrees from positive X-axis]')
linkaxes(axb,'x')

% Plot moving standard deviations
% figure
% axb(1) = tightSubplot(2,1,1);
% hold on
% for ii=1:length(iPix)
%     plot(t1,U1std(:,ii),'LineWidth',1.5,'Color',co0(ii,:));
% end
% for ii=1:length(iPix)
%     plot(V2.t,U2std(:,ii),'LineWidth',1.2,'Color',co(ii,:));
% end
% 
% 
% axb(2) = tightSubplot(2,1,2);
% hold on
% for ii=1:length(iPix)
%     plot(t1,U1std(:,ii),'LineWidth',1.5,'Color',co0(ii,:));
% end
% for ii=1:length(iPix)
%     plot(V2.t,U2std(:,ii),'LineWidth',1.2,'Color',co(ii,:));
% end

%% Plot 2nd order with various averaging lengths
meank = [ 3 7]; % 5 7];
Usmooth = zeros(length(V2.t),length(meank));
Theta_smooth = Usmooth;

Vpix = V2.Vz(60,180,:);
Upix = V2.Vx(60,180,:);
Uabs_raw = squeeze(sqrt(Vpix.^2 + Upix.^2));
Theta_raw = squeeze(atan2d(Vpix,Upix));
for jj = 1:length(meank)
    Usmooth(:,jj) = sqrt(movmean(Upix,meank(jj)).^2 + movmean(Vpix,meank(jj)).^2);
    Theta_smooth(:,jj) = atan2d(movmean(Vpix,meank(jj)),movmean(Upix,meank(jj)));
    
end

c2 = [0    0.5811    0.9633
         0    0.4470    0.7410
         0    0.2682    0.4446];
     
figure
axc(1)=tightSubplot(2,1,1);
plot(V2.t,Uabs_raw,'Color',[0.4 0.4 0.4],'LineWidth',2)
hold on
for jj = 1:length(meank)
    plot(V2.t,Usmooth(:,jj),'LineWidth',1.4,'Color',c2(jj,:))
end
axc(2)=tightSubplot(2,1,2);
plot(V2.t,Theta_raw,'Color',[0.4 0.4 0.4],'LineWidth',2)
hold on
for jj = 1:length(meank)
    plot(V2.t,Theta_smooth(:,jj),'LineWidth',1.4,'Color',c2(jj,:))
end
linkaxes(axc,'x')