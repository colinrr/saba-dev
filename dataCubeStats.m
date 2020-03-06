% Playing with a few ways to get temp statistics/visuals from the data cube

sDir        = '~/Kahuna/data/sabancaya_5_2018/image_exports/24A/spectral-calcs';
thermFile   = fullfile(sDir,'thermStats_2019-07-05_z131_x161_t598.mat');


 
load(thermFile)

%%
Tthresh = 250;
T = D.T;

% T = D.T-Tthresh;
% T(T<0)=0;
% T = T.^4;
% Get a horizontal "center of mass"
% T = (D.T - min(D.T(:))).^4;
% T = (D.T - Tthresh).^4;

% T = (D.T.^4 - min(D.T(:)).^4);

plumeC = squeeze(round( sum(reshape(1:numel(D.x),[1 numel(D.x) 1]).*T,2) ./ sum(T,2) ));

Tvar = squeeze(var(T,1,2))';

Lag0= C(1).Lag0;
Lmax = max(Lag0);
Lmin = min(Lag0);
i0 = Lmax-Lag0;
i1 = Lag0-Lmin;

nPts = numel(Tvar(1+i0(1):end-i1(1),1));
x = linspace(0,1,nPts)';
A2 = 1+i0 + x.*((nPts-i1) - (1+i0));
I = round(sub2ind(size(Tvar),A2,repmat(1:size(A2,2),[size(A2,1) 1])));

TvarA = Tvar(I);

Tmu = squeeze(mean(D.T,2))';

% Lag0= C(1).Lag0;
% Lmax = max(Lag0);
% Lmin = min(Lag0);
% i0 = Lmax-Lag0;
% i1 = Lag0-Lmin;

% nPts = numel(Tvar(1+i0(1):end-i1(1),1));
% x = linspace(0,1,nPts)';
% A2 = 1+i0 + x.*((nPts-i1) - (1+i0));
% I = round(sub2ind(size(Tvar),A2,repmat(1:size(A2,2),[size(A2,1) 1])));

TmuA = Tmu(I);
% tref = D.t(A2(:,refCol));

figure
subplot(2,1,1)
plot(D.z,TvarA(22,:))

subplot(2,1,2)
plot(D.z,TmuA(22,:))
%% 

checkZi = 86;

figure
subplot(2,1,1)
pcolor(D.t,D.z,plumeC)
shading flat

subplot(2,1,2)
pcolor(D.t,D.x,squeeze(D.T(checkZi,:,:)))
shading flat
hold on
plot(D.t,D.x(plumeC(checkZi,:)))