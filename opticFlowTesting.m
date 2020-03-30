clearvars -except D
close all

homeDir = '/Users/crrowell';
thermCube = fullfile(homeDir,'Kahuna/data/sabancaya_5_2018/image_exports/24A/thermCubeAnalysis/thermStats_2019-09-18_z641_x591_t1195.mat');

idx1 = 240;
dI   = 1;

method = 'classic+nl-fastp';
iSeq   = '1';

addpath(genpath(fullfile(homeDir,'Kahuna/data/sabancaya_5_2018/ijcv_flow_code/')));
%% Run and test
% cd ~/Kahuna/data/sabancaya_5_2018/ijcv_flow_code/

if ~exist('D','var')
    load(thermCube)
end

xcut = 1:250;
ycut = 1:250;
% Chop image size down around mask to run faster...
im1c = double(D.T(ycut,xcut,idx1));
im2c = double(D.T(ycut,xcut,idx1+dI));



myargs = {'lambda', 2};

tic
% uv = estimate_flow_demo('classic+nl-fastp', 4, 'middle-other', 'lambda', 3, 'pyramid_levels', 5);
uv1 = estimate_flow_interface(im1c,im2c,'classic+nl-fastp');
toc
tic
uv2 = estimate_flow_interface(im1c,im2c,'classic+nl-fastp', [], myargs); %, 'pyramid_levels', 5);
toc
% tic
% uv1 = estimate_flow_interface(im1,im2,'classic+nl-fastp');
% toc
% tic
% uv2 = estimate_flow_interface(im1,im2,'classic+nl-fastp', [], myargs); %, 'pyramid_levels', 5);
% toc


figure
tightSubplot(1,3,1)
plotThermVelocities(D.x(xcut),D.z(ycut),uv1(:,:,1),uv1(:,:,2),1)
tightSubplot(1,3,2)
plotThermVelocities(D.x(xcut),D.z(ycut),uv2(:,:,1),uv2(:,:,2),1)
tightSubplot(1,3,3)
% plotThermVelocities(D.x,D.z,uv(:,:,1),uv(:,:,2),1)
plotThermVelocities(D.x(xcut),D.z(ycut),uv2(:,:,1)-uv1(:,:,1),uv2(:,:,2)-uv1(:,:,2),.1)
