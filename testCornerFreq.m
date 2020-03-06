% Testing for slope change detection
clear all; %close all;

datadir   = '~/Kahuna/data/sabancaya_5_2018/image_exports/25/BI0525_big_explosion/';
thermdir  = fullfile(datadir,'interp-mat/');
specdata  = fullfile(datadir,'spectral-calcs/specT_1D_2019-0514_1633_w12_o6_n435.mat');

geom = fullfile(datadir,'mat/PTresults2/geometry.mat');
% Frame indices to use (row names in Tspectral table)
% idx    = [15 51 151 301 451 601 745]; 
idx    = {'451'}; 
winI   = 42;


%% Do the thing

load(specdata)
load(geom,'geom')


dat = table2struct(Tspectral(idx,:));
Ki  = dat.Wavenumber~=0;
K   = dat.Wavenumber(Ki);
K0  = zeros(size(dat.Pxx,2),1);
res = K0;

for ii = 1:size(dat.Pxx,2)
    Pxx = dat.Pxx(Ki,ii);
    Kint = linspace(min(log10(K)),max(log10(K)),numel(K));
    Pint = interp1(log10(K),log10(Pxx),Kint);

    % Find change points
    [ilt,res(ii)]=findchangepts(Pint,'MaxNumChanges',1,'Statistic','linear');
    K0(ii) = 10.^Kint(ilt);
end
%% 
wZ = dat.winZ-geom.Z0;

figure
loglog(dat.Wavenumber,dat.Pxx(:,winI))
hold on
plot(K0(winI)*[1 1],[min(dat.Pxx(winI,:)) max(dat.Pxx(winI,:))],'--k')

np=4;
figure
subplot(np,1,1)
plot(wZ,dat.Tmu)
axis tight
subplot(np,1,2)
plot(wZ,dat.specM)
axis tight
subplot(np,1,3)
plot(wZ,1./K0,'.-')
% hold on
axis tight
subplot(np,1,4)
% plot(wZ,dat.winD)
plot(wZ(1:end-1)+diff(wZ)/2,diff(dat.winD)./diff(wZ))
axis tight
grid on
