% clear all; close all
  
% datadir   = '/home/crowell/Kahuna/data/sabancaya_5_2018/';
% thermdir   = fullfile(datadir,'image_exports/24/1030_RT_explosion/');
% matDir  = fullfile(thermdir,'mat/');
% 
% params = fullfile(matDir,'params.mat');
% 
% gradDir = fullfile(matDir,'grad-scale/');
% parFile = fullfile(gradDir,'params.mat');
% 
% Idx  = [378:20:1300];
% rows = 6;
% cols = 8;
% bins = 190:2:424.74;
% %%
% binE = bins(1:end-1)+diff(bins);
% load(parFile)
% 
% hf = figure('position',[50 50 1300 800]);
% % ff = figure('position',[50 50 1300 1000]);
% 
% for ii=1:numel(Idx)
%     idx = Idx(ii);
%     load(fullfile(matDir,T.File{num2str(idx)}))
%     
%     % Histograms
%     N = histcounts(Frame(:),bins);
%     figure(hf)
%     tightSubplot(rows,cols,ii,0,0.02)
%     plot(binE,N,'.-')
%     set(gca,'YScale','log','YTickLabel',[])
%     axis tight
%     text(0.7,0.85,sprintf('%i',idx),'Units','normalized')
%     % Masks and frames
% end
%%
ITW = [378 420 700 900 1300 1567; ...
       260 255 250 240  230  230; ...
        10  10  10  10   10   10];
    
meth = 'pchip';
 
I = 378:1567;
T = interp1(ITW(1,:),ITW(2,:),I,meth);
dT = interp1(ITW(1,:),ITW(3,:),I,meth);

% figure
hold on
plot(ITW(1,:),ITW(2,:),'o')
% hold on
aa=plot(I,T);
plotLineError(I,T,dT,aa.Color)