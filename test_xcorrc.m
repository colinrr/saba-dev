% Test xcorrc a bit for vectorization
clear all; close all

dataDir = '~/Kahuna/data/sabancaya_5_2018/image_exports/24A/spectral-calcs/';
% thermFile   = fullfile(dataDir,'thermStats_2019-07-11_z151_x26_t598.mat');
thermFile   = fullfile(dataDir,'xcorrc_test_data.mat');

% Input for original xcorr
refCol = 1;
compCol = 2;
cutIdx  = [50 598-50];

xcflag = 'biased';
rep    = 1000;
perc   = 0.95;

% Filter
forder = 2;
flims  = [1/80];
ftype = 'high';

% Input for vectorized xcorr
compRows = [25 65];
%% Do the thing

load(thermFile)

% a = C.qDat(:,refCol);
% b = C.qDat(:,compCol);
% A = C.qDat;
% 
% 
% [y, bdu, bdl, auto_a, auto_b,lags,yy] = xcorrc(a,b,xcflag); %,rep,perc);
% 
% [y2, bdu2, bdl2, auto_a2, lags2,yy2] = xcorrcv(A,xcflag); %,rep,perc);
% 
% [y3,lags3] = xcorr(A,'coeff');

% [y, bdu, bdl, auto_a, auto_b,lags,yy] = xcorrc(a,xcflag);

% Synthetic
xs = linspace(0,2*pi,61);
ys = sin(xs);
ys1 = [zeros(150,1); ys'; zeros(150,1)];
ys2 = [zeros(185,1); ys'; zeros(115,1)];
ys3 = [zeros(194,1); ys'; zeros(106,1)];

a = ys1;
b = ys2;
A = [ys1 ys2 ys3];

tic
[y, bdu, bdl, auto_a, auto_b,lags,yy] = xcorrc(a,b,xcflag,rep,perc); %,rep,perc);
toc
tic
[y2, bdu2, bdl2, auto_a2, lags2,yy2] = xcorrcv(A,1,xcflag,rep,perc);
toc
[y3,lags3] = xcorr(A,xcflag);


% Preproc crap
% % Input for original xcorr
% a = detrend(D.Iint(refRow,cutIdx(1):cutIdx(2))','linear');
% B = detrend(D.Iint(compRows,cutIdx(1):cutIdx(2))','linear');
% t = D.t(cutIdx(1):cutIdx(2));
% 
% dt = mean(diff(t));
% fs = 1/dt;
% nyquist=fs/2;
% 
% % Filter
% [fb,fa]=butter(forder,flims./nyquist,ftype);
% a = filtfilt(fb,fa,a);
% B = filtfilt(fb,fa,B);
% b = B(:,1);
% 
% % Taper and pad?
% 
% [y, bdu, bdl, auto_a, auto_b,lags] = xcorrc(a,b,xcflag,rep,perc);


%% Plot

%     figure('position',[50 400 800 200])
%     plot(a)
%     hold on
%     plot(b)
%     axis tight
%     title('Input data')
% 
%     figure('position',[50 400 800 500])
%     subplot(2,1,1)
%     plot(lags,auto_a)
%     hold on
%     plot(lags,auto_b)
%     title('Original autocorr')
%     axis tight
% 
% %     figure('position',[50 400 800 200])
%     subplot(2,1,2)
%     plot(lags,y)
%     hold on
%     plot(lags,bdu)
%     plot(lags,bdl)
%     title('Original xcorrc and conf. intervals')
%     axis tight
    
figure
subplot(4,1,1)
plot(A,'LineWidth',3)
% hold on
% plot(b,'LineWidth',1.7)
title('Input Data')
axis tight

subplot(4,1,2)
plot(lags,y,'LineWidth',3)
hold on
plot(lags,yy,'LineWidth',1.7)
plot(lags,[bdu bdl])

title('xcorrc')
axis tight

subplot(4,1,3)
plot(lags2,y2(:,2),'LineWidth',3)
hold on
plot(lags2,yy2(:,2),'LineWidth',1.7)
title('xcorrcv')
plot(lags2,[bdu2(:,2) bdl2(:,2)])
axis tight

subplot(4,1,4)
plot(lags3,y3(:,2),'LineWidth',3)
axis tight
title('xcorr')
    
%     figure
%     plot(A)
%     hold on
%     plot(B)
%     plot(D.t,A(npad0+1:end-npad1))
%     hold on
%     plot(D.t,B(npad0+1:end-npad1))
