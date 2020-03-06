% Load and display particle size distribution info for Sabancaya
% clear all; close all

% ifile = '/Users/crrowell/Kahuna/data/sabancaya_5_2018/Result_C_Soriaux/Ash_sample_PSD.mat';
ifile = '~/Kahuna/data/sabancaya_5_2018/Result_C_Soriaux/Ash_sample_PSD.mat';



%% 

load(ifile)

figure
waterfall(centers_um,1:size(ash_samps_Mfrac,2),ash_samps_Mfrac')
colormap(flipud(winter))
xlabel('Size [{\mu}m]')
ylabel('Sample')
zlabel('Mass fraction')

hold on
plot3(centers_um,0*centers_um,mean(ash_samps_Mfrac,2)','k','LineWidth',2)
legend('Samples','Mean')
daspect([150 1 .1])
