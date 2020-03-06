% Get therm stats with height for each frame
% Equivelent to assuming steady plume, z=t

clearvars -except D
close all

dataDir   = '~/Kahuna/data/sabancaya_5_2018/image_exports/24A/';
thermFile = fullfile(dataDir,'thermCubeAnalysis/thermStats_2019-09-18_z641_x591_t1195.mat');

if ~exist('D','var')
    load(thermFile)
end

