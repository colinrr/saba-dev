% Rise diagrams and spectral analaysis
disp('')
disp('========= Spectral slope diagrams =========')
%% EVENT 24A
% fname  = 'specT_1D_2019-06-26_1213_w12_o8_n572.mat';
% datdir = '~/Kahuna/data/sabancaya_5_2018/image_exports/24/1030_RT_explosion/spectral-calcs';
% datfile = fullfile(datdir,fname);
specFile = fullfile(sDir,'specT_1D_2019-06-26_1213_w12_o8_n572.mat');

Tcax = [220 420];
% Scax = [-4 -1.8];
Scax = [-6 -0.7];
xl = [0 261.5];

pads = [0.15 0.1 0.1 0.05]; % row plots
pads = []; % column plots
% prfX = [20 45 90];
prfX = [30 60 90];
prfC = [1 1 1];

save_flag = true;
load_flag = false;

%% EVENT 25B
% fname   = ['RDspec_' datestr(now,'yyyy-mm-dd') '.mat'];
% fname  = 'RDspec_2019-03-14.mat';
% datdir = '~/Kahuna/data/sabancaya_5_2018/image_exports/25/BI0525_big_explosion/';
% matDir = fullfile(datDir,'mat/');
% datfile = fullfile(datdir,fname);
% specFile = fullfile(sDir,'specT_1D_2019-06-26_1641_w12_o8_n435.mat');
% 
% Tcax = [220 420];
% % Scax = [-4 -1.8];
% Scax = [-4 -2];
% xl = [16 103];
% 
% pads = [0.15 0.1 0.1 0.05];
% % prfX = [20 45 90];
% prfX = [30 60 90];
% prfC = [1 1 1];
% 
% save_flag = false;
% load_flag = false;
%% save data
% if save_flag % Run (or load) plumeSpec1D results and plotRiseDiagram
%     KS.t = Tspectral.VidTime;
%     KS.Z = Tspectral.winZ{end};
%     KS.S = Kolm;
%     
%     RD = D;
%     
%     save(datfile,'KS','RD','spec_params')
% end

%% load and plot data
load(geomf)
% load(datfile)

% Spectral data
load(specFile)
KS.t = Tspectral.VidTime;
KS.Z = Tspectral.winZ{end};
KS.S = Kolm;
%% Separate figures

figure('position',[20 20 800 450])
% axa=tightSubplot(1,1,1,0,0,pads);
pcolor(KS.t,KS.Z-geom.Z0,KS.S)
shading flat
view([0 0 1])
% colormap(parula(300))
colormap(CubeHelix(150))
caxis(Scax)
axis tight
cb2=colorbar('location','southoutside');
cb2.Label.String = 'Spectral slope';
% cb2.Position = [0.5600    0.1336    0.3750    0.0244];
% ylabel('Height [m]')
xlabel('time [s]')
ylabel('Height [m]')
% letterlabel3('b',[],14,'ilt');
hold on
xlim(xl)
zm=max(KS.S(:));
% plot3([prfX;prfX],[ylim' ylim' ylim'],zm*ones([2 3]),'--','LineWidth',2,'Color',prfC)
set(gca,'FontSize',12)

if save_flag
    ofig = sprintf('KolmDiagram_%s',datestr(now,'YYYY-dd-mm'));
%     printpdf(ofig,[25 15],figdir)
end

%% One figure, subplots as rows
% figure('position',[20 20 800 800])
% axa=tightSubplot(2,1,1,0,0,pads);
% set(gca,'FontSize',12)
% surf(RD.t,RD.z-RD.z0,RD.A)
% shading flat
% view([0 0 1])
% colormap(axa,thermgray(300))
% caxis(Tcax)
% axis tight
% cb = colorbar;%('location','west');
% cb.Label.String = 'T [K]';
% cb.Label.FontSize = 12;
% cb.Position([1 2 4]) = [0.92 0.55 0.375];
% hold on
% % yl=ylim';
% zm=max(RD.A(:));
% plot3([prfX;prfX],[ylim' ylim' ylim'],zm*ones([2 3]),'--','LineWidth',2,'Color',prfC)
% 
% ylabel('Height [m]')
% letterlabel3('a',[],14,'ilt','w');
% 
% axb=tightSubplot(2,1,2,0,0,pads);
% set(gca,'FontSize',12)
% surf(KS.t,KS.Z-RD.z0,KS.S)
% shading flat
% view([0 0 1])
% colormap(axb,gray(300))
% caxis(Scax)
% axis tight
% cb2=colorbar;%('location','west')
% cb2.Label.String = 'Spectral slope';
% cb2.Position([1 2 4]) = [0.92 0.125 0.375];
% ylabel('Height [m]')
% xlabel('Time [s]')
% letterlabel3('b',[],14,'ilt');
% hold on
% plot3([prfX;prfX],[ylim' ylim' ylim'],zm*ones([2 3]),'--','LineWidth',2,'Color',prfC)

%% Alternate (one figure, subplots as columns)

% figure('position',[20 20 1500 450])
% axa=tightSubplot(2,2,1,0,0,[],[],[3 1]);
% surf(RD.t,RD.z-RD.z0,RD.A)
% shading flat
% view([0 0 1])
% colormap(gca,thermgray(300))
% caxis(Tcax)
% axis tight
% xlabel('time [s]')
% cb = colorbar('location','southoutside');
% cb.Label.String = 'T [K]';
% cb.Label.FontSize = 12;
% cb.Position = [0.1300    0.1336    0.3750    0.0244];
% hold on
% % yl=ylim';
% xlim([16 103])
% zm=max(RD.A(:));
% plot3([prfX;prfX],[ylim' ylim' ylim'],zm*ones([2 3]),'--','LineWidth',2,'Color',prfC)
% 
% ylabel('Height [m]')
% set(gca,'FontSize',12)
% 
% % letterlabel3('a',[],14,'ilt','w');
% % figure('position',[20 20 800 450])
% axb=tightSubplot(2,2,2,0,0,[],[],[3 1]);
% surf(KS.t,KS.Z-RD.z0,KS.S)
% shading flat
% view([0 0 1])
% % colormap(parula(300))
% colormap(axb,CubeHelix(100))
% caxis(Scax)
% axis tight
% cb2=colorbar('location','southoutside');
% cb2.Label.String = 'Spectral slope';
% cb2.Position = [0.5600    0.1336    0.3750    0.0244];
% % ylabel('Height [m]')
% xlabel('time [s]')
% % letterlabel3('b',[],14,'ilt');
% hold on
% xlim([16 103])
% zm=max(KS.S(:));
% plot3([prfX;prfX],[ylim' ylim' ylim'],zm*ones([2 3]),'--','LineWidth',2,'Color',prfC)
% set(gca,'FontSize',12,'YTick',[])
% % 
